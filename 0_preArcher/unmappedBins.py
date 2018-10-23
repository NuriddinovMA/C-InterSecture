import sys
import timeit
from genome import Genome 

print 'Start genomes analysing...'
start_time = timeit.default_timer()
lines = sys.stdin.readlines()

for line in lines:
	parse = line.split()
	print '\tReading fasta...', parse[0]
	GNM = Genome(parse[0])
	i2l = GNM.idx2label
	for p in parse[1:]:
		r = int(p)
		print '\tSet resolution...', r
		GNM.setResolution(r)
		ubb = GNM.unmappedBasesBin
		f = open(parse[0] + '.%i.unmap' % r, 'w')
		for i in range(len(ubb)):
			for j in range(len(ubb[i])):
				print >> f, i2l[i],j*r,(j+1)*r-1,ubb[i][j]
		f.close()
		del ubb
		elp = timeit.default_timer() - start_time
		print '\t... end genome analysing for resolution %i, end time: %.2f sec' % (r,elp)
	del i2l
	elp = timeit.default_timer() - start_time
	print '\t... end genome analysing %.2f sec' % elp
elp = timeit.default_timer() - start_time
print '... full time %.2f sec' % elp