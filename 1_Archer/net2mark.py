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
for line in lines:
	l = line.strip()
	if l[0] != '#': pass
	else: file.append(l)
path = file[0]
for f in file[1].split():
	af.net2mark(path+f)
	elp = timeit.default_timer() - start_time
	print 'end converting:', file, elp
elp = timeit.default_timer() - start_time
print 'total end converting', elp