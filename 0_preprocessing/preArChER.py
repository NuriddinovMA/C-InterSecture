import sys
import timeit
import preArChER_func as prAf

Args = {
	'genome':'','chrom_order':'','sample_name':'','out_folder':'',
	'raw_contacts':'', 'norm_contacts':'', 'genome_bins':'',
	'statistic':'prc','resolution':10000, 'unmapped_bases':33,'coverage':1, 'max_distance':5
}

lines = sys.stdin.readlines()
for line in lines:
	parse = line.split('#')
	try:
		key = parse[0].strip().split()[0]
		args = parse[0].strip().split()[2].strip()
	except IndexError: pass
	else:
		print key, '=', args
		try: Args[key] = int(args)
		except ValueError: Args[key] = args
		except KeyError: pass


suffix = '/%s.%s.%ikb.%iN.%iC.%iMb' % (Args['sample_name'],Args['statistic'],Args['resolution']/1000,Args['unmapped_bases'],Args['coverage'],Args['max_distance'])
out_name = Args['out_folder'] + suffix
#
#This block indexes chromosomes and analyses bins by a unmapped_bases
#
print 'Step 0: data preparing...'
print '\tGenome analysis...'
start_time = timeit.default_timer()
l2i = prAf.ChromIndexing(Args['chrom_order'])
ubh = prAf.unmappedBasesBin(Args['genome'], Args['resolution'], l2i, Args['unmapped_bases'])
elp = timeit.default_timer() - start_time
print '\t... end genome analysing %.2f sec' % elp
#
#Reading bin order
#
print '\tBin labels reading...'
binIdxs = prAf.iBin2Label(Args['genome_bins'],l2i,Args['resolution'])
elp = timeit.default_timer() - start_time
print '\t...end bin labels reading %.2f sec' % elp
print '...end data preparing %.2f sec' % elp
#
#Reading raw contact matrix
#
print '\nStep 1: Raw matrix reading...'
print '\tDropped bins', len(ubh)
rawContactHash = prAf.iSparseMatrixReader(Args['raw_contacts'], binIdxs, ubh, coverage=Args['coverage'],raw=True)
print '\tDropped bins', len(ubh)
elp = timeit.default_timer() - start_time
print '\t%i contact analyzing for %.2f sec; memory sized: %.2f Mb' % (len(rawContactHash[1]), elp, 1.0*sys.getsizeof(rawContactHash)/1024/1024)

print '\nStep 2: Analyzing normalized marix for resolution %ikb... ' % (Args['resolution']/1000)
normContactHash = prAf.iSparseMatrixReader(Args['norm_contacts'], binIdxs, ubh, raw=False)
elp = timeit.default_timer() - start_time
print '\t%i normalized matrix reading tume: %.2f sec; memory sized: %.2f Mb' % (len(normContactHash),elp, sys.getsizeof(normContactHash )/1024.0/1024.0 )
#H1 = set( normContactHash.keys() )
#H2 = set( rawContactHash[0].keys() )
#H3 = set( rawContactHash[1].keys() )
#
#All contact hashed by distance
#
print '\nStep 3: Distance depended statistics... '
npContactBins = prAf.iContactBinStat(rawContactHash,normContactHash)
del rawContactHash
del normContactHash
elp = timeit.default_timer() - start_time
print '\t...genome bin statistic calculating end time  %.2f sec; sec memory sized: %.2f Mb' % ( elp, sys.getsizeof(npContactBins)/1024.0/1024.0 )
contactDistanceHash = prAf.iDistanceHash(npContactBins,Args['max_distance']*1e6/Args['resolution'])
del npContactBins
elp = timeit.default_timer() - start_time
print '...Analyzing end time %.2f sec; memory sized: %.2f Mb' % ( elp, sys.getsizeof(contactDistanceHash)/1024.0/1024.0 ) 
#
#Percentilyzing contact hashed by distance
#
print '\nStep 4: Contact transforming by %s statistic...' % Args['statistic']
totalContactList = prAf.iTotalContactListing(contactDistanceHash, Args['statistic'],Args['resolution'],out_name)
del contactDistanceHash
elp = timeit.default_timer() - start_time
print '... contact transforming end time %.2f sec; memory sized: %.2f Mb' % ( elp, sys.getsizeof(totalContactList)/1024.0/1024.0 )
#
#Writing contacts 
#
print '\nStep 5: Database writing...'
prAf.iBinContactWriter(out_name, totalContactList, l2i)
del totalContactList
del l2i
print '... writing contact end time %.2f sec' % elp
elp = timeit.default_timer() - start_time
print 'End: Full time for %ikb resolution %.2f sec' % (Args['resolution']/1000, elp)