# (c) 2012 Massachusetts Institute of Technology. All Rights Reserved
# Code written by: Anton Goloborodko (golobor@mit.edu),
# Maksim Imakaev (imakaev@mit.edu)



'''A Genome object contains the cached properties of a genome.

Glossary
--------

positional identifier : a number that can be decomposed into a tuple
    (chromosome index, base pair).

chromosome string label : the name of a chromosome.
    Examples: 'X', 'Y', '1', '22'.

chromosome index : a zero-based numeric index of a chromosome.
    For numbered chromosomes it is int(label) - 1, unless some of the
    chromosomes are absent. The chromosomes 'X', 'Y', 'M' are indexed
    after, the rest is indexed in alphabetical order.

concatenated genome : a genome with all chromosomes merged together
into one sequence.

binned genome : a genome splitted into bins `resolution` bp each.

binned concatenated genome : a genome with chromosomes binned and merged.
    GOTCHA: since the genome is binned FIRST and merged after that, the
    number of bins may be greater than (sum of lengths / resolution).
    The reason for this behavior is that the last bin in a chromosome
    is usually shorter than `resolution` but still counts as a full bin.
'''

from __future__ import absolute_import, division, print_function
import os
import glob
import re
import numpy
np = numpy
import warnings
import Bio.SeqIO
import Bio.SeqUtils
import Bio.Restriction
Bio.Restriction  # To shut up Eclipse warning
import joblib
#from scipy import weave
import logging
log = logging.getLogger(__name__)

h = logging.NullHandler()
log.addHandler(h)


class Genome(object):
    """A class to compute and cache various properties of genomic sequences."""

    def _memoize(self, func_name):
        '''Local version of joblib memoization.
        The key difference is that _memoize() takes into account only the
        relevant attributes of a Genome object (folder, name, gapFile,
        chrmFileTemplate) and ignores the rest.

        The drawback is that _memoize() doesn't check for the changes in the
        code of the function!'''
        if not hasattr(self, '_mymem'):
            self._mymem = joblib.Memory(cachedir=self.cacheDir)

        def run_func(readChrms, gapFile, chrmFileTemplate,
                     func_name, genomeName, *args, **kwargs):
            return getattr(self, func_name)(*args, **kwargs)

        mem_func = self._mymem.cache(run_func)

        def memoized_func(*args, **kwargs):
            return mem_func(
                self.readChrms, self.gapFile,
                self.chrmFileTemplate, func_name, self.folderName, *args, **kwargs)

        return memoized_func

    def setCentromeres(self,
            cntrStarts=None,
            cntrEnds=None,
            centromere_positions=None):
        """Set centromere positions.

        Parameters
        ----------

        cntrStarts, cntrEnds : np.arrays of int
            The arrays containing the starts and the ends of centromeres.
            Must contain chrmCount elements.

        centromere_positions : dict of (int, int) or np
            A dictionary with centromere positions.
            The keys are the chromosome string labels and the values are
            (centromereStart, centromereEnd). Can be supplied instead of
            `cntrStarts` and `cntrEnds`.
        """
        if ((centromere_positions is None)
              and not(cntrStarts is None)
              and not(cntrEnds is None)):

            cntrStarts = np.array(cntrStarts).astype(int)
            cntrEnds = np.array(cntrEnds).astype(int)
            assert(cntrStarts.size == self.chrmCount)
            assert(cntrEnds.size == self.chrmCount)
            self.cntrStarts = cntrStarts
            self.cntrEnds = cntrEnds

        elif (not(centromere_positions is None)
              and (cntrStarts is None)
              and (cntrEnds is None)):

            self.cntrStarts = numpy.zeros(self.chrmCount, int)
            self.cntrEnds = numpy.zeros(self.chrmCount, int)
            for label, (i, j) in list(centromere_positions.items()):
                chrm_idx = self.label2idx[label]
                self.cntrStarts[chrm_idx] = min(i, j)
                self.cntrEnds[chrm_idx] = max(i, j)
        else:
            raise Exception('Please provide either centromere_positions or ' +
                            'cntrStarts AND cntrEnds.')

        self.cntrMids = (self.cntrStarts + self.cntrEnds) // 2
        self.chrmArmLens = numpy.zeros(2 * self.chrmCount, int)
        self.chrmArmLens[0::2] = self.cntrMids
        self.chrmArmLens[1::2] = self.chrmLens - self.cntrMids
        lowarms = numpy.array(self.cntrStarts)
        higharms = numpy.array(self.chrmLens) - numpy.array(self.cntrEnds)
        self.maxChrmArm = max(lowarms.max(), higharms.max())

    def _parseGapFile(self):
        """Parse a .gap file to determine centromere positions.
        """
        gapPath = os.path.join(self.genomePath, self.gapFile)
        if not os.path.isfile(gapPath):
            log.warning(
                'Gap file not found! Set centromere positions to zero.\n')
            cntrStarts = numpy.zeros(self.chrmCount, int)
            cntrEnds = numpy.zeros(self.chrmCount, int)
        else:
            gapFile = open(gapPath).readlines()

            cntrStarts = -1 * numpy.ones(self.chrmCount, int)
            cntrEnds = -1 * numpy.zeros(self.chrmCount, int)

            for line in gapFile:
                splitline = line.split()
                if splitline[7] == 'centromere':
                    chrm_str = splitline[1][3:]
                    if chrm_str in self.label2idx:
                        chrm_idx = self.label2idx[chrm_str]
                        cntrStarts[chrm_idx] = int(splitline[2])
                        cntrEnds[chrm_idx] = int(splitline[3])
        self.setCentromeres(cntrStarts=cntrStarts, cntrEnds=cntrEnds)

    def createGapFile(self):
        """Create a gap file with the centromere positions.
        Use this method, if the genome you're using has no gap file.
        """
        gapPath = os.path.join(self.genomePath, self.gapFile)
        if os.path.isfile(gapPath):
            raise Exception('The gap file {0} already exists!'.format(gapPath))
        gapFile = open(os.path.join(self.genomePath, self.gapFile), 'w')
        for i in range(self.chrmCount):
            gapFile.write(
                '0\t{0}\t{1}\t{2}\t0\tN\t0\tcentromere\tno\n'.format(
                    (self.chrmFileTemplate % self.chrmLabels[i]).split('.')[0],
                    self.cntrStarts[i], self.cntrEnds[i]))

        gapFile.close()

    def _extractChrmLabel(self, fastaName):
        # First assume a whole filename as input (e.g. 'chr01.fa')
        _, fastaName = os.path.split(fastaName)
        regexp = self.chrmFileTemplate % ('(.*)')
        search_results = re.search(regexp, fastaName)
        # If not, assume that only the name is supplied as input (e.g. 'chr01')
        if search_results is None:
            regexp = self.chrmFileTemplate.split('.')[0] % ('(.*)')
            search_results = re.search(regexp, fastaName)

        if search_results is None:
            raise Exception(
                'The filename {} does not match the template {}.'.format(
                fastaName, self.chrmFileTemplate))

        chrm_label = search_results.group(1)

        # Remove leading zeroes.
        if chrm_label.isdigit():
            chrm_label = str(int(chrm_label))

        return chrm_label

    def _scanGenomeFolder(self):
        # Read chromosome IDs.
        self.chrmLabels = []
        filteredFastaNames = []
        for i in self.fastaNames:
            chrm = self._extractChrmLabel(i)
            if (not(self.readChrms)
                or
                (chrm.isdigit() and '#' in self.readChrms)
                or
                chrm in self.readChrms):

                self.chrmLabels.append(chrm)
                filteredFastaNames.append(i)
                log.debug('Convert %s FASTA filename to %s chromosome label, '
                              'store in the Genome object', i, chrm)
            else:
                log.debug('Convert %s FASTA filename to %s chromosome label, '
                              'discard', i, chrm)

        self.fastaNames = filteredFastaNames
        log.debug('The following FASTA files satisfy the readChrms variable '
                      '(={0}): {1}'.format(self.readChrms, self.fastaNames))

        if len(self.fastaNames) == 0:
            raise Exception('No Genome files at %s contain '
                            'the specified chromosomes' % self.genomePath)

        # Convert IDs to indices:
        # A. Convert numerical IDs.
        if not self.forceOrder:
            num_ids = [i for i in self.chrmLabels if i.isdigit()]
            log.debug('The chromosomes with numerical IDs: {0}'.format(num_ids))
            # Sort IDs naturally, i.e. place '2' before '10'.
            num_ids.sort(key=lambda x: int(re.findall(r'\d+$', x)[0]))

            self.chrmCount = len(num_ids)
            self.label2idx = dict(
                [(num_ids[i], int(i)) for i in range(len(num_ids))])
            self.idx2label = dict(
                [(int(i), num_ids[i]) for i in range(len(num_ids))])
        else:
            self.chrmCount = 0
            self.label2idx = {}
            self.idx2label = {}

        # B. Convert non-numerical IDs. Give the priority to XYM over the rest.
        if self.forceOrder:
            nonnum_ids = self.readChrms
        else:
            nonnum_ids = sorted([i for i in self.chrmLabels if not i.isdigit()])

        log.debug('The chromosomes with non-numerical IDs: {0}'.format(
            nonnum_ids))
        if not self.forceOrder:
            for i in ['M', 'Y', 'X']:
                if i in nonnum_ids:
                    nonnum_ids.pop(nonnum_ids.index(i))
                    nonnum_ids.insert(0, i)
        for i in nonnum_ids:
            self.label2idx[i] = self.chrmCount
            self.idx2label[self.chrmCount] = i
            self.chrmCount += 1

        # Sort fastaNames and self.chrmLabels according to the indices:
        self.chrmLabels = list(zip(*sorted(list(self.idx2label.items()),
                                      key=lambda x: x[0])))[1]
        self.fastaNames.sort(
            key=lambda path: self.label2idx[self._extractChrmLabel(path)])
        log.debug('The genome folder is scanned successfully.')


    def __init__(self,
                 genomePath,
                 gapFile='gap.txt',
                 chrmFileTemplate='chr%s.fa',
                 readChrms=['#', 'X', 'Y', 'M'],
                 cacheDir='default', forceOrder=False):
        '''
        A class that stores cached properties of a genome. To initialize,
        a Genome object needs FASTA files with chromosome sequences.
        For the basic definitions please refer to the module documentation.

        Parameters
        ----------
        (for the constructor method)

        genomePath : str
            The path to the folder with the FASTA files.

        gapFile : str
            The path to the gap file relative to genomePath.

        chrmFileTemplate : str
            The template of the FASTA file names.

        readChrms : list of str
            The list with the string labels of chromosomes to read from the
            genome folder. '#' stands for chromosomes with numerical labels
            (e.g. 1-22 for human). If readChrms is empty then read all
            chromosomes.

        cacheDir : directory, or "default" for caching in the genomePath

        Attributes
        ----------

        chrmCount : int
            the total number of chromosomes.

        chrmLabels : list of str
            a list of chromosomal IDs sorted in ascending index order.

        fastaNames : list of str
            FASTA files for sorted in ascending index order of respective
            chromosomes.

        genomePath : str
            The path to the folder with the genome.

        folderName : str
            The string identifier of the genome, the name of the last folder in
            the path.

        label2idx : dict
            a dictionary for conversion between string chromosome labels and
            zero-based indices.

        idx2label : dict
            a dictionary for conversion between zero-based indices and
            string chromosome labels.

        seqs : list of str
            a list of chromosome sequences. Loads on demand.

        chrmLens : list of int
            The lengths of chromosomes.

        chrmArmLens : list of int
            The lengths of chromosomal arms including the centromeric regions.

        maxChrmLen : int
            The length of the longest chromosome.

        cntrStarts : array of int
            The start positions of the centromeres.

        cntrMids : array of int
            The middle positions of the centromeres.

        cntrEnds : array of int
            The end positions of the centromeres.

        The following attributes are calculated after setResolution() is called:

        resolution : int
            The size of a bin for the binned values.

        chrmLensBin : array of int
            The lengths of chromosomes in bins.

        chrmArmLensBin : array of int
            The lengths of chromosomal arms in bins.

        chrmStartsBinCont : array of int
            The positions of the first bins of the chromosomes in the
            concatenated genome.

        chrmEndsBinCont : array of int
            The positions of the last plus one bins of the chromosomes in the
            concatenated genome.

        chrmBordersBinCont : array of int
            The positions of the chromosome borders in the concatenated genome.
            chrmBordersBinCont equals to chrmStartsBinCont with the total
            number of bins + 1 appended.

        chrmIdxBinCont : array of int
            The index of a chromosome in each bin of the concatenated genome.

        posBinCont : array of int
            The index of the first base pair in a bin in the concatenated
            genome.

        cntrMidsBinCont : array of int
            The position of the middle bin of a centromere in the concatenated
            genome.

        chrmArmBordersBinCont: array of int
            The position of the chromosome arm borders in the concatenated
            genome.

        GCBin : list of arrays of float
            % of GC content of bins in individual chromosomes.

        unmappedBasesBin : list of arrays of int
            Number of bases with N's for each bin

        mappedBasesBin : list of arrays of int
            Number of sequenced bases for each bin

        binSizesbp : list of arrays of int
            Size of each bin. Is less than *resolution* for the last bin only.

        The following attributes are calculated after setEnzyme() is called:

        enzymeName : str
            The restriction enzyme used to find the restriction sites.

        rsites : list of arrays of int
            The indices of the first base pairs of restriction fragments
            in individual chromosomes.

        rfragMids : list of arrays of int
            The indices of the middle base pairs of restriction fragments
            in individual chromosomes.

        rsiteIds : array of int
            The position identifiers of the first base pairs of restriction
            fragments.

        rsiteMidIds : array of int
            The position identifiers of the middle base pairs of restriction
            fragments.

        rsiteChrms : array of int
            The indices of chromosomes for restriction sites in corresponding
            positions of rsiteIds and rsiteMidIds.
        '''
        # Set the main attributes of the class.
        self.genomePath = os.path.abspath(os.path.expanduser(genomePath))
        self.forceOrder = forceOrder
        self.singleFile = False
        if os.path.isfile(self.genomePath):
            self.singleFile = self.genomePath
            self.genomePath = os.path.dirname(self.genomePath)

        if cacheDir == "default":
            cacheDir = self.genomePath
        self.cacheDir = cacheDir
        self.folderName = os.path.split(self.genomePath)[-1]
        self.readChrms = list(readChrms)
        self.gapFile = gapFile
        self.chrmFileTemplate = chrmFileTemplate

        log.debug('Initialize a Genome object genomePath=%s, readChrms=%s, '
                    'gapFile=%s, chrmFileTemplate=%s', self.genomePath,
                    self.readChrms, self.gapFile, self.chrmFileTemplate)

        # Scan the folder and obtain the list of chromosomes.
        if not self.singleFile:
            self.fastaNames = [os.path.join(self.genomePath, i)
                               for i in glob.glob(os.path.join(
                                                  self.genomePath,
                                                  self.chrmFileTemplate % ('*',)))]
            log.debug('Scan genome folder: {0}'.format(self.genomePath))
            log.debug('FASTA files are found: {0}'.format(self.fastaNames))
            self._scanGenomeFolder()
        else:
            import pyfaidx
            with open(self.singleFile, 'r') as f:
                self.fastaNames = list(rec.name for rec in pyfaidx.Fasta(self.singleFile))
            self.chrmLabels = self.fastaNames
            self.chrmCount = len(self.chrmLabels)
            self.label2idx = dict([(self.chrmLabels[i], i) for i in range(self.chrmCount)])
            self.idx2label = dict([(i, self.chrmLabels[i]) for i in range(self.chrmCount)])

        # Get the lengths of the chromosomes.
        self.chrmLens = self.getChrmLen()
        self.maxChrmLen = max(self.chrmLens)
        # FragIDmult is used in (chrm, frag) -> fragID conversion.
        self.fragIDmult = self.maxChrmLen + 1000

        # Parse a gap file and mark the centromere positions.
        self._parseGapFile()

    def getChrmLen(self):
        # At the first call redirects itself to a memoized private function.
        self.getChrmLen = self._memoize('_getChrmLen')
        return self.getChrmLen()

    def _getChrmLen(self):
        return numpy.array([len(self.seqs[i])
                            for i in range(0, self.chrmCount)])

    def getGCBin(self, resolution):
        # At the first call the function rewrites itself with a memoized
        # private function.
        self.getGCBin = self._memoize('_getGCBin')
        return self.getGCBin(resolution)

    def _getGCBin(self, resolution):
        GCBin = []
        for chrm in range(self.chrmCount):
            chrmSizeBin = int(self.chrmLens[chrm] // resolution) + 1
            GCBin.append(numpy.ones(chrmSizeBin, dtype=numpy.float))
            for j in range(chrmSizeBin):
                GCBin[chrm][j] = self.getGC(
                    chrm, j * int(resolution), (j + 1) * int(resolution))
        return GCBin

    def getUnmappedBasesBin(self, resolution):
        # At the first call the function rewrites itself with a memoized
        # private function.
        self.getUnmappedBasesBin = self._memoize('_getUnmappedBasesBin')
        return self.getUnmappedBasesBin(resolution)

    def _getUnmappedBasesBin(self, resolution):
        unmappedBasesBin = []
        for chrm in range(self.chrmCount):
            chrmSizeBin = int(self.chrmLens[chrm] // resolution) + 1
            unmappedBasesBin.append(numpy.ones(chrmSizeBin, dtype=numpy.float))
            for j in range(chrmSizeBin):
                unmappedBasesBin[chrm][j] = self.getUnmappedBases(
                    chrm, j * int(resolution), (j + 1) * int(resolution))
        return unmappedBasesBin

    def getRsites(self, enzymeName):
        # At the first call redirects itself to a memoized private function.
        self.getRsites = self._memoize('_getRsites')
        return self.getRsites(enzymeName)

    def _getRsites(self, enzymeName):
        '''Returns: tuple(rsites, rfrags)
        Finds restriction sites and mids of rfrags for a given enzyme
        '''

        # Memorized function
        rsites = []
        rfragMids = []

        if type(enzymeName) == int:
            for i in range(self.chrmCount):
                rsites.append(np.r_[np.arange(0, self.chrmLens[i]-100, enzymeName), self.chrmLens[i]])
        else:
            enzymeSearchFunc = eval('Bio.Restriction.%s.search' % enzymeName)
            for i in range(self.chrmCount):
                rsites.append(numpy.r_[
                    0, numpy.array(enzymeSearchFunc(self.seqs[i].seq)) + 1,
                    len(self.seqs[i].seq)])

        for i in range(self.chrmCount):
            rfragMids.append((rsites[i][:-1] + rsites[i][1:]) // 2)

        # Remove the first trivial restriction site (0)
        # to equalize the number of points in rsites and rfragMids.
        for i in range(len(rsites)):
            rsites[i] = rsites[i][1:]

        return rsites, rfragMids

    @property
    def seqs(self):
        if not hasattr(self, "_seqs"):
            if not self.singleFile:
                self._seqs = []
                for i in range(self.chrmCount):
                    self._seqs.append(Bio.SeqIO.read(open(self.fastaNames[i]),
                                                     'fasta'))
            else:
                import pyfaidx
                fasta_records = pyfaidx.Fasta(self.singleFile)
                self._seqs = []
                for name in self.chrmLabels:
                    # monkey patch to behave like SeqIO object for now
                    seq = Bio.Seq.Seq(str(fasta_records[name]))
                    seq.seq = seq
                    self._seqs.append(seq)
        return self._seqs

    def setResolution(self, resolution):
        """Set the resolution of genome binning and calculate the following
        attributes:

        resolution : int
            The size of a bin for the binned values.

        chrmLensBin : array of int
            The lengths of chromosomes in bins.

        chrmArmLensBin : array of int
            The lengths of chromosomal arms in bins.

        chrmStartsBinCont : array of int
            The positions of the first bins of the chromosomes in the
            concatenated genome.

        chrmEndsBinCont : array of int
            The positions of the last plus one bins of the chromosomes in the
            concatenated genome.

        chrmBordersBinCont : array of int
            The positions of the chromosome borders in the concatenated genome.
            chrmBordersBinCont equals to chrmStartsBinCont with the total
            number of bins + 1 appended.

        chrmIdxBinCont : array of int
            The index of a chromosome in each bin of the concatenated genome.

        posBinCont : array of int
            The index of the first base pair in a bin in the concatenated
            genome.

        cntrMidsBinCont : array of int
            The position of the middle bin of a centromere in the concatenated
            genome.

        chrmArmBordersBinCont: array of int
            The position of the chromosome arm borders in the concatenated
            genome.

        GCBin : list of arrays of float
            % of GC content of bins in individual chromosomes.

        unmappedBasesBin : list of arrays of int
            Number of bases with N's for each bin

        mappedBasesBin : list of arrays of int
            Number of sequenced bases for each bin

        binSizesbp : list of arrays of int
            Size of each bin. Is less than *resolution* for the last bin only.
        """

        if (resolution == -1) and hasattr(self, "resolution"):
            for i in ["chrmLensBin", "chrmArmLensBin",
                      "chrmBordersBinCont",
                      "chrmStartsBinCont",
                      "chrmEndsBinCont", "numBins",
                      "chrmIdxBinCont", "posBinCont",
                      "cntrMidsBinCont", "chrmArmBordersBinCont",
                      "GCBin",
                      "unmappedBasesBin", "binSizesBp",
                      "mappedBasesBin", "resolution"]:
                exec("del self.%s" % i)
            return

        self.resolution = int(resolution)

        # Bin chromosomes.
        self.chrmLensBin = self.chrmLens // self.resolution + 1
        self.chrmBordersBinCont = numpy.r_[0, numpy.cumsum(self.chrmLensBin)]
        self.chrmStartsBinCont = numpy.r_[0, numpy.cumsum(
            self.chrmLensBin)[:-1]]
        self.chrmEndsBinCont = numpy.cumsum(self.chrmLensBin)
        self.numBins = self.chrmEndsBinCont[-1]

        self.chrmIdxBinCont = numpy.zeros(self.numBins, int)
        self.binStartPosCont = numpy.zeros(self.numBins, int)
        self.binEndPosCont = numpy.zeros(self.numBins, int)
        self.posBinCont = numpy.zeros(self.numBins, int)
        
        for i in range(self.chrmCount):
            st = self.chrmStartsBinCont[i]
            end = self.chrmEndsBinCont[i]
            self.chrmIdxBinCont[st:end] = i
            self.binStartPosCont[st:end] = resolution * np.arange(end - st, dtype = int)
            self.binEndPosCont[st:end-1] = self.binStartPosCont[st+1:end]
            self.binEndPosCont[end-1] = self.chrmLens[i]

        
        for i in range(self.chrmCount):
            self.posBinCont[
                self.chrmStartsBinCont[i]:self.chrmEndsBinCont[i]] = (
                    self.resolution
                    * numpy.arange(-self.chrmStartsBinCont[i]
                                   + self.chrmEndsBinCont[i]))

        # Bin centromeres.
        self.cntrMidsBinCont = (self.chrmStartsBinCont
                                + self.cntrMids // self.resolution)
        self.chrmArmBordersBinCont = numpy.zeros(
            self.chrmCount * 2 + 1, dtype=numpy.int)
        self.chrmArmBordersBinCont[1::2] = self.cntrMidsBinCont
        self.chrmArmBordersBinCont[2::2] = self.chrmEndsBinCont
        self.chrmArmLensBin = self.chrmArmBordersBinCont[1:] - \
            self.chrmArmBordersBinCont[:-1]

        # Bin GC content.
        self.GCBin = self.getGCBin(self.resolution)
        self.unmappedBasesBin = self.getUnmappedBasesBin(self.resolution)
        self.binSizesBp = []
        for i in range(self.chrmCount):
            chromLen = self.chrmLens[i]
            cur = [self.resolution for _ in range(chromLen // self.resolution)]
            cur.append(chromLen % self.resolution)
            self.binSizesBp.append(numpy.array(cur))
        self.mappedBasesBin = [numpy.array(i[0] * (100. - i[1]) / 100, int)
                for i in zip(self.binSizesBp, self.unmappedBasesBin)]

    def getUnmappedBases(self, chrmIdx, start, end):
        "Calculate the percentage of unmapped base pairs in a region."
        seq = self.seqs[chrmIdx][start:end]
        if hasattr(seq, "seq"):
            myseq = seq.seq
        else:
            myseq = seq 
        if len(myseq) == 0:
            return 0.0
        else:
            return (100.0 * (myseq.count('N') + myseq.count('n'))
                    / float(len(myseq)))

    def getGC(self, chrmIdx, start, end):
        """Calculate the GC content of the mapped part of a region. If there
        are no mapped base pairs, return 50%.
        """
        seq = self.seqs[chrmIdx][start:end]
        if hasattr(seq, "seq"):
            myseq = seq.seq
        else:
            myseq = seq 
        overall_GC = Bio.SeqUtils.GC(myseq)
        unmapped_content = self.getUnmappedBases(chrmIdx, start, end)

        if unmapped_content == 100.0:
            return -1.0
        else:
            corrected_GC = overall_GC * 100.0 / (100.0 - unmapped_content)
            return corrected_GC

    def clearCache(self):
        '''Delete the cached data in the genome folder.'''
        if hasattr(self, '_mymem'):
            self._mymem.clear()

    def setEnzyme(self, enzymeName, rsites=None):
        """Apply a specified restriction enzyme to the genomic sequences and
        calculate the positions of restriction sites.

        Parameters
        ----------

        enzymeName : str
            The name of the restriction enzyme from Biopython Restriction module.

        rsites : list of np.arrays
            Use only if your restriction enzyme is not present in Biopython.
            Contains positions of restriction sites for each chromosome.
            By default is None and ignored.

        The following attributes are set with this method:
        --------------------------------------------------

        enzymeName : str
            The restriction enzyme used to find the restriction sites.

        rsites : list of arrays of int
            The indices of the first base pairs of restriction fragments
            in individual chromosomes.

        rfragMids : list of arrays of int
            The indices of the middle base pairs of restriction fragments
            in individual chromosomes.

        rfragLens : list of arrays of int
            The lengths of restriction fragments in bp.

        rsiteIds : array of int
            The position identifiers of the first base pairs of restriction
            fragments.

        rsiteMidIds : array of int
            The position identifiers of the middle base pairs of restriction
            fragments.

        rsiteChrms : array of int
            The indices of chromosomes for restriction sites in corresponding
            positions of rsiteIds and rsiteMidIds.

        chrmStartsRfragCont : array of int
            The absolute indices of first restriction fragments in each
            chromosome.

        chrmEndsRfragCont : array of int
            The absolute indices of last restriction fragments in each
            chromosome.

        cntrMidsRfrag : array of int
            The relative indices of restriction fragments containing the
            centromere midpoint.

        cntrMidsRfragCont : array of int
            The absolute indices of restriction fragments containing the
            centromere midpoint.

        chrmBordersRfragCont : array of int
            The indices of restriction fragments delimiting each chromosome.
        """

        if not (rsites is None):
            if enzymeName in Bio.Restriction.__dict__:
                raise Exception('Since you are using a custom restriction '
                    'enzyme, you cannot use a name of an existing enzyme.')
            if len(rsites) != self.chrmCount:
                raise Exception('Please, specify restriction sites for each '
                                'chromosome.')

            rsites = [np.sort(i) for i in rsites]
            rfragMids = []
            for i in range(self.chrmCount):
                if rsites[i][0] != 0:
                    rsites[i] = np.hstack([[0], rsites[i]])
                if rsites[i][-1] != 0:
                    rsites[i] = np.hstack([rsites[i], len(self.seqs[i].seq)])
                rfragMids.append((rsites[i][:-1] + rsites[i][1:]) / 2)
                rsites[i] = rsites[i][1:]
            self.rsites = rsites
            self.rfragMids = rfragMids
        else:
            self.rsites, self.rfragMids = self.getRsites(enzymeName)

        self.enzymeName = enzymeName
        self.rfragLens = [
            numpy.diff(numpy.r_[0, i]) for i in self.rsites]
        self.chrmEndsRfragCont = numpy.cumsum([len(i) for i in self.rsites])
        self.chrmBordersRfragCont = numpy.r_[0, self.chrmEndsRfragCont]
        self.chrmStartsRfragCont = self.chrmBordersRfragCont[:-1]
        self.cntrMidsRfrag = [
            numpy.searchsorted(self.rsites[i], self.cntrMids[i]) - 1
            for i in range(self.chrmCount)]
        self.cntrMidsRfragCont = self.cntrMidsRfrag + self.chrmStartsRfragCont

        self.numRfrags = self.chrmEndsRfragCont[-1]

        self.rsiteIds = numpy.concatenate(
            [self.rsites[chrm] + chrm * self.fragIDmult
             for chrm in range(self.chrmCount)])

        self.rfragMidIds = numpy.concatenate(
            [self.rfragMids[chrm] + chrm * self.fragIDmult
             for chrm in range(self.chrmCount)])

        self.rsiteChrms = numpy.concatenate(
            [numpy.ones(len(self.rsites[chrm]), int) * chrm
             for chrm in range(self.chrmCount)])

        assert (len(self.rsiteIds) == len(self.rfragMidIds))

    def hasEnzyme(self):
        return hasattr(self, "enzymeName")

    def splitByChrms(self, inArray):
        return [inArray[self.chrmStartsBinCont[i]:self.chrmEndsBinCont[i]]
                for i in range(self.chrmCount)]

    def upgradeMatrix(self, oldGenome):
        """Checks if old genome can be upgraded to new genome by truncation.
        If not, returns an array that can be used
        to upgrade chromosome positions.
        If upgrade not possible, raises an exception.

        Parameters
        ----------
        old Genome : Genome, or label2idx dictionary
            old genome from which upgrade is done

        Returns
        -------
        None : upgrade is possible by truncating chromosomes >= chromNum
        upgradeIndex : ndarray  upgrade is possible
        by newChrom = upgradeMatrix[oldChrom]

        Raises an exception when upgrade is not possible
        """

        if isinstance(oldGenome, Genome):
            oldGenome = oldGenome.idx2label
        if True in [i not in list(oldGenome.values())
                    for i in list(self.idx2label.values())]:
            difference = [i for i in list(self.idx2label.values(
                )) if i not in list(oldGenome.values())]
            raise Exception("Genome upgrade is not possible: " +
                            repr(difference) + " are chromosomes"
                            " that are missing in the old genome")
        if False not in [oldGenome[i] == self.idx2label[i]
                         for i in list(self.idx2label.keys())]:
            return None
        oldLabelToIdx = dict([(oldGenome[i], i) for i in list(oldGenome.keys())])
        convertingArray = numpy.zeros(len(list(oldGenome.keys())), dtype=int) - 1
        for i in list(self.idx2label.values()):
            convertingArray[oldLabelToIdx[i]] = self.label2idx[i]
        return convertingArray

    def checkReadConsistency(self, chromosomes, positions):
        """

        """
        chromSet = set(chromosomes)
        if 0 not in chromSet:
            warnings.warn("Chromosome zero not found! Are you using"
                          " zero-based chromosomes?", UserWarning)
        if max(chromSet) >= self.chrmCount:
            raise Exception("Chromosome number %d exceeds expected"
                                " chromosome count %d" %
                                (max(chromSet), self.chrmCount))
        if max(chromSet) < self.chrmCount - 1:
            warnings.warn("More chromosomes in the genome (%d)  than we got"
                          " (%d) ! Are you using proper genome?" %
                          (self.chrmCount, max(chromSet) - 1))
        maxpositions = self.chrmLens[chromosomes]
        check = positions > maxpositions
        if check.any():  # found positions that exceeds chromosme length
            inds = numpy.nonzero(check)[0]
            inds = inds[::len(inds) // 10]
            for i in inds:
                raise Exception("Position %d on chrm %d exceeds "
                                    "maximum positions %d" % (
                        chromosomes[i], positions[i],
                        self.chrmLens[chromosomes[i]])
                                    )

    def getRfragAbsIdxs(self, rfragIds):
        """Convert restriction fragment IDs to absolute fragment indices.

        Parameters
        ----------
        rfragIds: array of int
            IDs of fragments, calculated as
            fragIDmult * chromosome + fragment midpoint

        Returns
        -------
        rfragAbsIdxs: array of int
            absolute indices of restriction fragments in a concatenated
            genome.
        """

        if not self.hasEnzyme():
            raise Exception('Please set the restriction enzyme first.')

        rfragAbsIdxs = numpy.searchsorted(self.rfragMidIds, rfragIds)
        # map the fragment absolute indices back to the relative indices and
        # check every 100th for correctness
        assert (rfragIds[::100] - self.rfragMidIds[rfragAbsIdxs[::100]]).sum() == 0

        return rfragAbsIdxs

    def getFragmentDistance(self, fragments1, fragments2, enzymeName):
        """returns distance between fragments
        measured in... fragments. (neighbors = 1, etc. )"""
        if not hasattr(self, "rfragMidIds"):
            self.setEnzyme(enzymeName)
        frag1ind = numpy.searchsorted(self.rfragMidIds, fragments1)
        frag2ind = numpy.searchsorted(self.rfragMidIds, fragments2)
        distance = numpy.abs(frag1ind - frag2ind)
        del frag1ind, frag2ind
        ch1 = fragments1 // self.fragIDmult
        ch2 = fragments2 // self.fragIDmult
        distance[ch1 != ch2] = 1000000
        return distance

    def getPairsLessThanDistance(self, fragments1, fragments2,
                                 cutoffDistance, enzymeName):
        from . import numutils
        """returns all possible pairs (fragment1,fragment2)
        with fragment distance less-or-equal than cutoff"""
        if not hasattr(self, "rfragMidIds"):
            self.setEnzyme(enzymeName)
        f1ID = numpy.searchsorted(self.rfragMidIds, fragments1)
        f2ID = numpy.searchsorted(self.rfragMidIds, fragments2)

        assert (fragments1[::100] - self.rfragMidIds[f1ID[::100]]).sum() == 0
        assert (fragments2[::100] - self.rfragMidIds[f2ID[::100]]).sum() == 0

        fragment2Candidates = numpy.concatenate(
            [f1ID + i for i in (list(range(-cutoffDistance, 0)) +
                                list(range(1, cutoffDistance + 1)))])
        fragment1Candidates = numpy.concatenate(
            [f1ID for i in (list(range(-cutoffDistance, 0)) +
                            list(range(1, cutoffDistance + 1)))])
        mask = numutils.arrayInArray(fragment2Candidates, f2ID)

        fragment2Real = fragment2Candidates[mask]
        fragment1Real = fragment1Candidates[mask]
        return  (self.rfragMidIds[fragment1Real],
                 self.rfragMidIds[fragment2Real])

    def _parseFixedStepWigAtSmallResolution(self, filename, resolution):
        """Internal method for parsing fixedStep wig file
        and averaging it over every kb"""
        if resolution > 10000:
            raise ValueError("Please use parseWigFile instead for resolutions >10 kb")
        myfilename = filename
        if os.path.exists(filename) == False:
            raise Exception("File not found!")
        M = self.maxChrmLen
        Mkb = int(M // resolution + 1)
        chromCount = self.chrmCount
        data = numpy.zeros(Mkb * self.chrmCount, np.double)
        resolution = int(resolution)
        if "X" in self.chrmLabels:
            useX = True
            Xnum = self.label2idx["X"] + 1  # wig uses zero-based counting
        else:
            useX = False
            Xnum = 0

        if "Y" in self.chrmLabels:
            useY = True
            Ynum = self.label2idx["Y"] + 1
        else:
            useY = False
            Ynum = 0

        if "M" in self.chrmLabels:
            useM = True
            Mnum = self.label2idx["M"] + 1
        else:
            useM = False
            Mnum = 0

        from .fastExtensions.fastExtensionspy import readWigFile  # @UnresolvedImport
        readWigFile(myfilename, data, chromCount,
                            useX, useY, useM,
                            Xnum, Ynum, Mnum, Mkb, resolution)

        datas = [data[i * Mkb:(i + 1) * Mkb] for i in range(self.chrmCount)]
        for chrom, track in enumerate(datas):
            if track[self.chrmLens[chrom] // resolution + 2:].sum() != 0:
                raise Exception("Genome mismatch: entrees "
                                    "in wig file after chromosome end!")
        datas = [numpy.array(i[:self.chrmLens[chrom] // resolution +
            1]) for chrom, i in enumerate(datas)]
        return datas

    def parseFixedStepWigAtSmallResolution(self, filename, resolution=5000):
        "Returns averages of a fixedStepWigFile for all chromosomes"
        # At the first call the function rewrites itself with a memoized
        # private function.
        self.parseFixedStepWigAtSmallResolution = self._memoize(
            '_parseFixedStepWigAtSmallResolution')
        return self.parseFixedStepWigAtSmallResolution(filename,
                                                    resolution=resolution)

    def _parseBigWigFile(self, filename, resolution=1000,
                         divideByValidCounts=False):

        if resolution > 10000:
            raise ValueError("Please use parseWigFile instead for resolutions >10 kb")
        import bx.bbi.bigwig_file
        from bx.bbi.bigwig_file import BigWigFile  # @UnresolvedImport

        """
        Internal method for parsing bigWig files, updated
        """
        data = []
        if type(filename) == str:
            bwFile = BigWigFile(open(filename))
        else:
            bwFile = BigWigFile(filename)
        print("parsingBigWigFile", end=' ')
        assert isinstance(bwFile, bx.bbi.bigwig_file.BigWigFile)  # @UndefinedVariable

        for i in range(self.chrmCount):
            chrId = "chr%s" % self.idx2label[i]
            print(chrId, end=' ')
            totalCount = int(numpy.ceil(self.chrmLens[i] / float(resolution)))
            values = numpy.zeros(totalCount, float)
            step = 500
            for i in range(totalCount // step):
                beg = step * i
                end = min(step * (i + 1), totalCount * resolution)
                summary = bwFile.summarize(chrId, beg *
                    resolution, end * resolution, end - beg)
                if summary is None:
                    continue
                stepValues = summary.sum_data
                stepCounts = summary.valid_count
                if divideByValidCounts == True:
                    stepValues = stepValues / stepCounts
                    stepValues[stepCounts == 0] = 0
                values[beg:end] = stepValues
            if values.sum() == 0:
                raise  Exception("Chromosome {0} is absent in bigWig"
                                     " file!".format(chrId))
            data.append(values)

        return data

    def parseBigWigFileAtSmallResolution(self, filename, resolution=1000,
                        divideByValidCounts=False):
        """
        Parses bigWig file using bxPython build-in method "summary".
        Does it by averaging values over "resolution" long windows.

        If window has less than lowCountCutoff valid valies, it is discarded

        Parameters
        ----------

        filename : str or file object
            Incoming bigWig file
        lowCountCutoff : int, < resolution
            Ignore bins with less than cutoff valid bases
        resolution : int
            Find average signal over these bins
        divideByValidCounts : bool
            Divide  total coverage of the kb bin.

        Returns
        -------
        List of numpy.arrays with average values for each chromosomes
        Length of each array is ceil(chromLens / resolution)
        """

        # At the first call the function rewrites itself with a memoized
        # private function.
        self.parseBigWigFile = self._memoize('_parseBigWigFile')
        return self.parseBigWigFile(filename,
                                    resolution, divideByValidCounts)
    def parseAnyWigFile(self, filenames, control=None,
                    wigFileType="Auto", functionToAverage=np.log, internalResolution=1000):
        """
        1. Calculates total value of a track over internalResolution (1kb) sized subbins

        2. If multiple files are given, aggregates them at this level by taking a sum

        3. Aggregated reads are divided by control within each subbin

        4. We take functionToAverage (by default - log) of all non-zero subbins in the bin

        5. We takes an average of all non-zero subbins within a bin

        6. We divide each track by the median of the track

        Parameters
        ----------

        filenames : str or list of str
            list of filenames to aggregate data from
        control : str or None
            Control track
        wigFileType : "bigwig", "wig" or "auto"
            Is suppled file bigwig, or fixed step wig
        functionToAverage : function
            Averages f(1kb values in a bin) over current resolution bin
        internalResolution : int, divider of self.resolution, multiple of 1kb
            Internal resolution to summarize a bigwig or wig file over

        Returns
        -------

        List of by chromosome averages of a track at self.resolution resolution



        Only fixed-step wig files and bigWig files are supported!!!
        Import from fixedStep wig files is very fast,
        however is not the most reliable.

        for VariableStep files use wigToBigWig utility to convert
        them to bigWig format first.
        To use it you will also need to have fetchChromSizes script.

        Then you just run
        $bash fetchChromSizes hg18 > hg18.chrom.sizes
        $./wigToBigWig myWig.wig hg18.chrom.sizes myWig.bigWig

        And you enjoy your favourite bigWig.

        BigWig import is implemented using bx-python module.
        It is normally very fast; however, it has a bug at low resolutions.
        I have an ugly workaround for it (chopping the quiery
        into many smaller pieces), but I hope
        that they actually fix this bug.

        Anyways, check their repo on BitBucket, maybe they've
        fixed my issue # 39 :)
        https://bitbucket.org/james_taylor/bx-python/overview
        """

        if type(filenames) == str:
            filenames = [filenames]
        filenames = [os.path.abspath(i) for i in filenames]
        wigFileType = wigFileType.lower()

        def loadFile(name, wigFileType=wigFileType):
            """Choosing right method to load wigfile"""

            if os.path.exists(name) == False:
                raise IOError("\n Wig file not found : %s " %
                    (os.path.abspath(name)))

            if wigFileType == "auto":
                ext = os.path.splitext(name)[1]
                if ext == "":
                    raise Exception("Wig file has no extension. \
                    Please specify it's type")
                elif ext.lower() == ".wig":
                    op = open(name)
                    if "fi" not in [op.readline()[:2] for _ in range(5)]:
                        raise Exception("Cannot read non variable-step wig \
                        files! Please use wigToBigWig utility. See docstring \
                        of this method.")
                    wigFileType = "wig"
                elif ext.lower() == ".bigwig":
                    wigFileType = "bigwig"
                else:
                    raise Exception("Unknown extension of wig file: %s" %
                        ext)

            if wigFileType == "wig":
                data = self.parseFixedStepWigAtSmallResolution(
                    name, resolution=internalResolution)
            elif wigFileType == "bigwig":
                data = self.parseBigWigFileAtSmallResolution(name, resolution=internalResolution,
                    divideByValidCounts=True)
            else:
                raise Exception("Wrong type of wig file : %s" %
                    wigFileType)
            return data

        "1. Calculates total value of a track over internalResolution (5kb) sized subbins"
        data = [loadFile(i) for i in filenames]

        "2. If multiple files are given, aggregates them at this level by taking a sum"
        for otherdata in data[1:]:
            for i in range(len(otherdata)):
                data[0][i] += otherdata[i]
        data = data[0]

        "3. loading control"
        if control is not None:
            controlData = loadFile(control)

        if self.resolution % internalResolution != 0:
            raise Exception("Cannot parse wig file at resolution \
            that is not a multiply of internal resolution ({0})".format(internalResolution))

        resultByChromosome = []
        for chrom, value in enumerate(data):
            value = np.array(value)
            if control is not None:
                chromControl = np.asarray(controlData[chrom])
                vmask = value != 0
                cmask = chromControl != 0
                keepmask = vmask * cmask  # keeping only bins with non-zero reads in data/control
                vmasksum, cmasksum = vmask.sum(), cmask.sum()
                # Comparing number of non-zero bins in a control and non-control group
                if max(vmasksum, cmasksum) / (1. * min(vmasksum, cmasksum)) \
                > 1.3:
                    warnings.warn("\nBig deviation: number of non-zero \
                    data points: %s, control points:%s."
                    % (vmasksum, cmasksum))
                value[-keepmask] = 0
                "3. Aggregated reads are divided by control within each subbin"
                value[keepmask] = value[keepmask] / chromControl[keepmask]

            # Making a linear array into a matrix (rows - bins, columns - subbins within a bin)
            # possibly appending some extra zeros at the end
            value.resize(self.chrmLensBin[chrom] * (
                self.resolution / internalResolution))

            value.shape = (-1, self.resolution / internalResolution)
            if value.mean() == 0:
                raise Exception("Chromosome {0} contains zero data in wig \
                file(s) {1}".format(self.idx2label[chrom], filenames))

            "4. We take functionToAverage (by default - log) of all non-zero subbins in the bin"
            mask = value == 0
            value[-mask] = functionToAverage(value[-mask])

            valuesum = np.sum(value, axis=1)
            masksum = np.sum(mask == False, axis=1)
            "5. We takes an average of all non-zero subbins within a bin"
            valuesum[masksum == 0] = 0
            vmask = valuesum != 0
            valuesum[vmask] /= masksum[vmask]
            valuesum[-vmask] = 0
            # setting all unknown points to the median of known points
            resultByChromosome.append(valuesum)
        return resultByChromosome