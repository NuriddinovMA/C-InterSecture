import sys
import math
import timeit
import ArcherStatFunction as FS
import matplotlib
import matplotlib.pyplot as plt
#import matplotlib.colors as colors
import numpy as np

reload(FS)
import DataDifs
#reload(DataDifs)
from DataDifs import *

start_time = timeit.default_timer()
p = -1, 101
p1up = 25,101
p1dw = 0,75
p2up = 0,75
p2dw = 25,101
dp = 33
vl = 0
s = 125,125
if method == 'abs': T = True
else: T = False
xy = [ [],[] ]
xyr = [ [ ([],[]) for j in range(5) ], [ ([],[]) for j in range(5) ] ]
for f in range(len(allCon)):
	for sp in range(1):#range(len(allCon[f])):#
		print name[f][sp]
		print 'start metric calc'
		M = FS.Metric(allCon[f][sp],8,T)
		for i in range(4):
			fn = open(path_out+'%s_%s.%i.metric.bedGraph' % (file_path[f][:-1],name[f][sp],i),'w')
			for m in M: print >> fn, '%s\t%i\t%i\t%f' % (m[0],m[1]*Resolution,m[1]*Resolution+Resolution,m[2+i])
			fn.close()
		del M
		elp = timeit.default_timer() - start_time
		print '...all calculating total time:', elp
	#del allCon[f]