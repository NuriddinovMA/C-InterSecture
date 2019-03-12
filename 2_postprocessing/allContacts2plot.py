import os
import sys
import timeit
import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import post_func as psf
reload(psf)

print 'Step 1: Initialization'
start_time = timeit.default_timer()

Args = {
	'contact_path':'','samples':[],'contact_files':[],
	'frame':[], 'use_loci':False, 'loci': [], 'statistics': False,
	'chrom_path':'','chrom_sizes':{}, 
	'use_synblocks':False, 'synblocks_path':'', 'synblocks_files':[],
	'out_path':'','out_names':[]
}

lines = sys.stdin.readlines()
for line in lines:
	parse = line.split('#')
	try:
		key = parse[0].strip().split()[0]
		args = parse[0].strip().split()[2:]
	except IndexError: pass
	else:
		if key == 'frame': Args[key] = [ int(s) for s in args ]
		elif key == 'contact_files': Args[key].append( [ s for s in args ] )
		elif key == 'samples' or key == 'out_names' or key == 'loci': Args[key].extend( [ s for s in args ] )
		elif key == 'chrom_sizes':
			if Args.has_key(key) == True: Args[key].update( dict([ s.split(':') for s in args ]))
			else: Args[key] = dict([ s.split(':') for s in args ])
		elif key == 'synblocks_files': Args[key].append( dict([ s.split(':') for s in args ]) )
		elif key == 'use_synblocks' or key == 'statistics' or key == 'use_pre' or key == 'use_loci': Args[key] = psf.boolean(args[0])
		else:
			try: Args[key] = args[0]
			except KeyError: pass

if Args['use_synblocks'] == True and len(Args['synblocks_files']) != len(Args['samples']): 
	print "Error! Number of 'synblock_files' pairs don't equal to number of samples!"
	exit()

if len(Args['contact_files']) != len(Args['samples']): 
	print "Error! Number of 'contact_files' strings don't equal to number of samples!"
	exit()

if len(Args['out_names']) == 0: Args['out_names'] = Args['samples']
elif len(Args['out_names']) == len(Args['samples']): pass
elif len(Args['out_names']) == 1: 
	Args['out_names'] *= len(Args['samples'])
	print "Attention! All output files share the prefix!"
else:
	print "Error! Number of 'out_names' must correspond to number of samples"

L = []

if Args['use_loci'] == True:
	for l in Args['loci']:
		parse = l.split(':')
		key = parse[0]
		try:
			chrm = parse[1]
			start,end = parse[2].split('-')
			start,end = int(start),int(end)
			L.append((key,chrm,start,end) )
		except IndexError: pass
	Args['loci'] = {}
	for l in L:
		if Args['loci'].has_key(l[0]) == True: Args['loci'][l[0]].append( l[1:] )
		else: Args['loci'][l[0]] = [ l[1:], ]
else: Args['loci'] = []

del L
for key in Args.keys(): print '\t%s =' % key, Args[key]

c = 0
for sm in range(len(Args['samples'])):
	if Args['contact_files'][sm][0] == 'all' or Args['contact_files'][sm][0][:3] == 'key': files = os.listdir( '%s/%s/' % (Args['contact_path'],Args['samples'][sm]))
	else: files = [ a for a in Args['contact_files'][sm] ]
	if Args['contact_files'][sm][0][:3] == 'key':
		for fi in range(len(files)-1,-1,-1):
			if files[fi].find(Args['contact_files'][sm][0][4:]) == -1: del files[fi]
	for file in files:
		s = file.split('.')
		fname = '%s/%s/%s' % (Args['contact_path'],Args['samples'][sm],file)
		if Args['chrom_sizes'].has_key(s[0]) == True and s[1] == 'prc' and s[-1] == 'liftContacts':
			elp = timeit.default_timer() - start_time
			print '\tstart contact analizying %s, %s, %.2f' % (Args['samples'][sm], file, elp)
			resolution = int(s[2][:-2])*1000
			order_path = '%s/%s' % (Args['chrom_path'],Args['chrom_sizes'][s[0]])
			Order = psf.ChromIndexing(order_path)
			loci = Args['loci'][s[0]]
			allCon = psf.readContacts(fname,Order,resolution)
			lociCon = psf.filterLoci(allCon,resolution,loci)
			del allCon
			elp = timeit.default_timer() - start_time
			print '\tend contact reading: %.2f' % elp

			for l in range(len(loci)):
				elp = timeit.default_timer() - start_time
				print '\tstart draw loci %s:%i-%i, time: %.2f' % (loci[l][0],loci[l][1],loci[l][2],elp)
				Map = psf.drawMap(lociCon[l],loci[l],resolution,1,1,3)
				Ln = len(Map[6])
				fig, ax = plt.subplots()
				ax1 = plt.subplot2grid((11,11),(0,0),colspan=8,rowspan=8)
				sc1 = plt.scatter(Map[1], Map[0], s=Map[4], c=Map[5], cmap = 'coolwarm', vmin=-100, vmax=100)
				sc2 = plt.scatter(Map[0], Map[1], s=Map[2], c=Map[3], cmap = 'autumn_r', vmin=0, vmax=100)
				ax1.set_aspect('equal',anchor='W')
				step = int(pow(Ln,0.5))
				MarkLines = [ i for i in range(0,Ln,step) ]
				for i in MarkLines:
					plt.axhline(y=i, linewidth=2, color='g', ls='--',alpha=0.5)
					plt.axvline(x=i, linewidth=2, color='g', ls='--',alpha=0.5)
					plt.text(i, i, i+1, withdash=True, size = 12)
					plt.text(i, i, i+1, withdash=True, size = 12)
					plt.text(-3, i,i+1, withdash=True, size = 12)
					plt.text(i, -2, i+1, withdash=True, size = 12)
					plt.text(Ln+4, i, i+1, withdash=True, size = 24)
					plt.text(i, Ln+6, i+1, withdash=True, size = 24)
				if Ln < 30: Ln = 30
				size = Ln/144.0
				plt.xticks( MarkLines, Map[6][::step], rotation=90, size = 24 )
				plt.yticks( MarkLines, Map[6][::step], size = 24 )
				ax1.invert_yaxis()
				ax1.xaxis.tick_top()
				
				del Map
				ax2 = plt.subplot2grid((12,12),(0,9),rowspan=8)
				ticks = [j for j in range(-100,105,10)]
				cb = plt.colorbar(sc1, cax=ax2, ticks = ticks)#, pad=0.05, fraction=0.05)
				cb.ax.tick_params(labelsize=24, length = 10, width = 2)
				cb.ax.set_ylabel('observed - control percentile score',fontsize=24 )
				ax2.set_aspect('auto')
				
				ax4 = plt.subplot2grid((12,12),(9,0),colspan=8)
				ticks = [j for j in range(0,105,5)]
				cb = plt.colorbar(sc2, cax=ax4, ticks = ticks, orientation='horizontal')#, pad=0.1, fraction=0.05)
				cb.ax.tick_params(labelsize=24, length = 10, width = 2)
				cb.ax.set_xlabel('observed percentile score',fontsize=24 )
				ax4.set_aspect('auto')
				
				Map = psf.drawStat(1,1,3)
				ax3 = plt.subplot2grid((12,12),(0,11),rowspan=8)
				x = np.arange(101)
				y = np.zeros(101)
				sc3 = plt.scatter(y[::2]+.5, x[::2], s=Map[1][::2], c=y[::2]-25, cmap = 'coolwarm', vmin=-100, vmax=100)
				plt.xlim([0,1])
				plt.xticks([.5,],['Liftovered',], size = 24)
				plt.yticks(x[::5],x[::5], size = 24)
				plt.text(1.5,63,'percentile deviation', size=24, rotation = 90)
				#ax4.invert_xaxis()
				ax3.yaxis.tick_right()
				ax3.set_aspect('auto')
				
				ax5 = plt.subplot2grid((12,12),(11,0),colspan=8)
				x = np.arange(101)
				y = np.zeros(101)
				sc5 = plt.scatter(x[::2], y[::2]+.5, s=Map[0][::2], c=y[::2]+50, cmap = 'autumn', vmin=0, vmax=100)
				plt.ylim([0,1])
				plt.yticks([.5,],['Observed',], size = 24)
				plt.xticks(x[::5],x[::5], size = 24)
				plt.xlabel('percentile deviation', size=24)
				ax5.set_aspect('auto')
				
				map = sc1.get_figure()
				map.set_size_inches((48*size,48*size))
				#plt.tight_layout()
				plt.savefig( Args['out_path'] + '/%s_%s.%s.%s.%i-%i.map.png' % (Args['out_names'][sm],s[0],s[1],loci[l][0],loci[l][1],loci[l][2]),dpi=100)
				map.clear()
				elp = timeit.default_timer() - start_time
				print '\t\tend plotting %.2f' % elp
elp = timeit.default_timer() - start_time
print '...end contact comparing', elp