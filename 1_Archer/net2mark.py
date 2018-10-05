import sys
import timeit
import ArChER_func as af

reload(FS)

#
#This utility transforms net-file to the more usefull (for ArChER) mark-file 
#
lines = sys.stdin.readlines()
path = lines[0].split('=')[1].strip()
file = []
for i in range (1,len(lines)):
	file.append( lines[i].split('=')[1].strip() )
print 'Step 1: net file reading...'
start_time = timeit.default_timer()
for i in file:
	af.net2mark(path+i) 
elp = timeit.default_timer() - start_time
print 'end net file converting', elp

