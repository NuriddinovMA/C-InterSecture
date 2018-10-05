import sys
import timeit
from genome import Genome 
import numpy as np
import preArChER_func as pf

reload(pf)

lines = sys.stdin.readlines()
genome_path = lines[0].split('=')[1].strip()
genome_file = lines[1].split('=')[1].strip()
contact_path = lines[2].split('=')[1].strip()
contact_file = lines[3].split('=')[1].strip()
Resolution = int(lines[4].split('=')[1].strip())
polyN = int(lines[5].split('=')[1].strip())
Cover_thold = int(lines[6].split('=')[1].strip())
Limit = int(lines[7].split('=')[1].strip())*1e6/Resolution
stat = lines[8].split('=')[1].strip()
out_path = lines[9].split('=')[1].strip()
out_name = lines[10].split('=')[1].strip()

#
#This block indexes chromosomes and analyses bins by a mappability
#
print 'Step 0: chromosome preparing...'
print '\tReading fasta'
start_time = timeit.default_timer()
GNM = Genome(genome_path+genome_file)
elp = timeit.default_timer() - start_time
print '\tend file reading', elp
l2i = GNM.label2idx
fname = out_path + out_name + '.%iN.%iC.%ikb.%iMb' % (polyN,Cover_thold,Resolution/1000,Limit*Resolution/1000000)
print '\tWrite chromosome index to file', fname
File = open(fname, 'w')
for key in l2i: print >> File, '%s\t%i' % (key, l2i[key])
File.close()
l2i.update( GNM.idx2label )
print '\t... end chromosome indexing', elp
GNM.setResolution(Resolution)
print '\tGenome N-trackt analysis...'
ubb = GNM.unmappedBasesBin
ubh = {}
for i in range(len(ubb)):
	for j in range(len(ubb[i])):
		if ubb[i][j] > polyN: ubh[(i,j)] = int(ubb[i][j])
del ubb
elp = timeit.default_timer() - start_time
print '\t... end genome analysing', elp

#
#Reading bin order
#
print 'Step 1: contact matrix read...', elp
File = open(contact_path + contact_file + '/raw/%i/' % Resolution + contact_file+'_%i_abs.bed' % Resolution,'r')
lines = File.readlines()
File.close()

iRB2L = pf.iReadBin2Label(lines,l2i, Resolution)
del lines
#
#Reading raw contact matrix
#
File = open(contact_path + contact_file +'/raw/%i/' % Resolution + contact_file +'_%i.matrix' % Resolution,'r')
lines = File.readlines()
File.close()
elp = timeit.default_timer() - start_time
print '\tend file reading', elp, 'memory sized:', sys.getsizeof(lines)
print '... contact reading end time', elp
#
#Filtring bins by coverage and read number
#
print 'Step 2: contact analyzing... %ikb... ' % (Resolution/1000)
flt_Bin_hash = pf.iRawMatrixFilter(lines, iRB2L, Cover_thold, ubh)
print len(ubh),
ubh.update(flt_Bin_hash)
print len(ubh)
Bin_hash = pf.iRawMatrixRead(lines, iRB2L, ubh)
del lines
del flt_Bin_hash
del ubh
elp = timeit.default_timer() - start_time
print "\t%i contact analyzing:" % Bin_hash[2], elp, 'memory sized:', sys.getsizeof(Bin_hash)
#
#Reading normalized contact matrix
#
print 'Step 3: Reading normalized marix for resolution %ikb... ' % (Resolution/1000)
File = open(contact_path + contact_file +'/iced/%i/' % Resolution + contact_file+'_%i_iced.matrix' % Resolution,'r')
lines = File.readlines()
File.close()
IMR = pf.IcedMatrixRead(lines, iRB2L)
H1 = set( IMR.keys() )
H2 = set( Bin_hash[0].keys() )
H3 = set( Bin_hash[1].keys() )
#print len(H1), len(H2), len(H3)
#print len(H1 | H2),len(H1 | H3),len(H2 | H3)
#print len(H1 & H2),len(H1 & H3),len(H2 & H3)
#print len(H1 - H2),len(H1 - H3),len(H2 - H3)
#print len(H2 - H1),len(H3 - H1),len(H3 - H2)
del lines
#
#All contact hashed by distance
#
Locus_contact_np = np.array(pf.iContactBinStat(Bin_hash,IMR,Limit), dtype=np.float32)
elp = timeit.default_timer() - start_time
print "\t...genome binnig end time", elp, 'memory sized:', sys.getsizeof(Locus_contact_np)
Con_VC = pf.iDistanceHash(Locus_contact_np)
elp = timeit.default_timer() - start_time
print "...Normalization end time", elp, 'memory sized:', sys.getsizeof(Con_VC)
#
#Percentilyzing contact hashed by distance
#
print "Step 4:  contact..."
del Locus_contact_np
Con_prc_total = []
File = open(out_path + out_name+'.%iN.%iC.%ikb.%iMb' % (polyN,Cover_thold,Resolution/1000,Limit*Resolution/1000000),'w')
Keys = Con_VC.keys()
Keys.sort()
for key in Keys:
	if key == -100: print >> File, 'Interchromosome\t%i' % len(Con_VC[key])
	else: print >> File, '%ikb\t%i' % (key*Resolution/1000, len(Con_VC[key])) 
	if stat == 'abs':
		Con_prc_list = pf.PercentileStatistics(Con_VC[key])
		Con_prc_total += pf.Percentilize(Con_VC[key], Con_prc_list)
	else:
		Con_prc_total += pf.iAbsoluteContacts(Con_VC[key])
	Con_VC[key] = 0
	elp = timeit.default_timer() - start_time
	print "\tpercentile calculated for %ikb distance" % (key*Resolution/1000), elp, 'memory sized:', sys.getsizeof(Con_prc_total)
del Con_VC
del Keys
File.close()
elp = timeit.default_timer() - start_time
print "... percentilize end time", elp, 'memory sized:', sys.getsizeof(Con_prc_total)
#
#Writing contacts 
#
print "Step 5: Database writing..."
pf.BinContactWriting(out_path+out_name+'.%iN.%iC.%ikb.%iMb' % (polyN,Cover_thold,Resolution/1000,Limit*Resolution/1000000), Con_prc_total, l2i)
print "... writing contact end time", elp
elp = timeit.default_timer() - start_time
print "End: Full time for %ikb" % (Resolution/1000), elp