import timeit
import genome as gnm
import preArChER_func as pf

reload(pf)

lines = sys.stdin.readlines()
path = lines[0].split('=')[1].strip()
file = lines[1].split('=')[1].strip().split()
suffix = '.fa'
print 'Start chromosome indexing...'
start_time = timeit.default_timer()
File = open (path+file+suffix,'r')
lines = File.readlines()
File.close()
elp = timeit.default_timer() - start_time
print '\tend file reading', elp
print '\tReads analysis...'
Chr_Ind = pf.ChrIndexing(lines)
pf.PrintChrIndex(Chr_Ind, path+file)
elp = timeit.default_timer() - start_time
print '\t... end chromosome indexing', elp
