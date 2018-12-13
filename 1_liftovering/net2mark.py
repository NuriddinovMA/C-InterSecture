import sys
import timeit
import ArChER_func as af

#reload(af)

#
#This utility transforms net-file to the more usefull (for ArChER) mark-file 
#
print 'Step 1: net file reading...'
start_time = timeit.default_timer()
file = []
lines = sys.stdin.readlines()
for line in lines:
	l = line.strip()
	if l[0] != '#': file.append(l)
	else: pass
path = file[0]
for flist in file[1:]:
	for f in flist.split():
		nP = af.netParser(path+f)
		n2p = af.net2pre(nP, path+f)
		del nP
		af.pre2mark(n2p, path+f)
		del n2p
		elp = timeit.default_timer() - start_time
		print 'end converting:', file, elp
elp = timeit.default_timer() - start_time
print 'total end converting', elp