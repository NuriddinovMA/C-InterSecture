import sys
import math
import timeit
import matplotlib
import matplotlib.pyplot as plt
#import matplotlib.colors as colors
import numpy as np
import scipy.stats as sc
import post_func as psf

reload(psf)

start_time = timeit.default_timer()

resolution = 50000
base_dir = '//fishminja/share/ArChER/'
file_path = 'allC/FbHM/','allC/FbHG/'
chrom_path = 'chrSizes/hg38.chr.sizes', 'chrSizes/hg38.chr.sizes', 'chrSizes/mm10.chr.sizes'#'chrSizes/galGal5.chr.sizes'
suf = '.prc.%ikb.33N.1C.5Mb.balanced.allContacts' % (resolution/1000)
sp_name = ['IMR90','IMR90']#
ln = len(sp_name)
repeat = 3

a = (0,100)
sp = 'FbHM_','FbHG_'
type = 'pp'
limit = '%i-%i_' % a
r = [a[0]-100,a[1]]
bin = 50#r[1] - r[0]
s = 1

for i in range(2):#ln):
	x = []
	y = []
	xr = []
	yr = []
	fname = base_dir+file_path[i]+sp_name[i]+suf 
	orname = base_dir+chrom_path[i]
	elp = timeit.default_timer() - start_time
	print '\tstart %s contact read total time: %.2f sec' % (sp_name[i], elp)
	Order = psf.ChromIndexing(orname)
	allCon = psf.readContacts(fname,Order,resolution,short=True)
	elp = timeit.default_timer() - start_time
	print '\t...end contact reading time: %.2f' % elp
	print '\tstart contact stat time'
	for k in allCon:
		if 0 < k[2] <= 20:
			x.append(k[0])
			y.append(k[1])
	elp = timeit.default_timer() - start_time
	print '\t...end contact stat time: %.2f' % elp
	print '\tstart randome contact stat'
	for j in range(repeat):
		elp = timeit.default_timer() - start_time
		print '\t\tstart %i repeat contact randomizing: %.2f sec' % (j, elp)
		randCon = psf.randomizeContacts(allCon)
		elp = timeit.default_timer() - start_time
		print '\t\t...end contact randomizing: %.2f sec' % elp
		print '\t\tstart random contact stat time'
		for k in range(len(randCon)-1,-1,-1):
			if 0 < randCon[k][2] <= 20:
				xr.append(randCon[k][0])
				yr.append(randCon[k][1])
			del randCon[k]
		del randCon
		elp = timeit.default_timer() - start_time
		print '\t...end random contact stat time: %.2f' % elp
	del allCon
	elp = timeit.default_timer() - start_time
	print '\t...total end contact stat time: %.2f' % elp
	plt.cla()
	H = [[],[]]
	#print bin
	#print 'spearman:'
	#print sc.spearmanr(x,y)
	print 'hist'
	H[0] = np.histogram2d(x,y,bins = (bin,bin), range=[(0,100),(0,100)], normed=True)
	del x
	del y
	H[1] = np.histogram2d(xr,yr,bins = (bin,bin), range=[(0,100),(0,100)], normed=True)
	del xr
	del yr
	M = np.log2(np.transpose(H[0][0]/H[1][0]))#
	M1 = [M[j] for j in range(bin-1,-1,-1)]
	del M
	im = plt.imshow(M1,cmap='viridis',vmin=0,vmax=5)
	plt.colorbar(im)
	plt.xticks( range(0,bin+1,10), [H[0][1][k] for k in range(0,bin+1,10)], fontsize = 10, rotation=45 )
	plt.yticks( range(0,bin+1,10), [H[0][2][k] for k in range(bin,-1,-10)], fontsize = 10 )
	#plt.show()
	plt.savefig( base_dir + sp[i] + type + '.png',dpi=400)
	plt.clf()

