import os
import sys
import timeit
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import post_func as psf
reload(psf)

print 'Step 1: Initialization'
start_time = timeit.default_timer()

Args = {
	'contact_path':'','samples':[],'contact_files':[],
	'frame':[], 'use_loci': False,'loci': [], 'statistics': False, 'metric': 'pbad',
	'chrom_path':'','chrom_sizes':{}, 
	'use_synblocks':False, 'synblocks_path':'', 'synblocks_files':[],
	'out_path':'','out_names':''
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
		elif key == 'contact_files' or key == 'loci': Args[key].append( [ s for s in args ] )
		elif key == 'samples' or key == 'out_names': Args[key] = [ s for s in args ]
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
elif len(Args['out_names']) == 1: print "Attention! All output files share the prefix!"
else:
	print "Error! Number of 'out_names' must correspond to number of samples"

if Args['use_loci'] == True:
	L = []
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
	del L
else: Args['loci'] = []
for key in Args.keys(): print '\t%s =' % key, Args[key]

elp = timeit.default_timer() - start_time
c = 125
if Args['metric'] == 'pearsone' or  Args['metric'] == 'spearman': tail = 'equ'
else: tail = 'prc'

elp = timeit.default_timer() - start_time
print 'Step 2: Analyzing', elp
for sm in range(len(Args['samples'])):
	if Args['contact_files'][sm][0] == 'all' or Args['contact_files'][sm][0][:3] == 'key': files = os.listdir( '%s/%s/' % (Args['contact_path'],Args['samples'][sm]))
	else: files = [ a for a in Args['contact_files'][sm] ]
	if Args['contact_files'][sm][0][:3] == 'key':
		for fi in range(len(files)-1,-1,-1):
			if files[fi].find(Args['contact_files'][sm][0][4:]) == -1: del files[fi]
	colorList = psf.colorList(len(files))
	for file in files:
		s = file.split('.')
		fname = '%s/%s/%s' % (Args['contact_path'],Args['samples'][sm],file)
		if Args['chrom_sizes'].has_key(s[0]) == True and s[-1] == 'liftContacts' and s[1] == tail:
			elp = timeit.default_timer() - start_time
			print '\tstart contact analizying %s, %s, %.2f' % (Args['samples'][sm], file, elp)
			resolution = int(s[2][:-2])*1000
			order_path = '%s/%s' % (Args['chrom_path'],Args['chrom_sizes'][s[0]])
			Order = psf.ChromIndexing(order_path)
			allCon = psf.readContacts(fname,Order,resolution)
			elp = timeit.default_timer() - start_time
			print '\tend contact reading: %.2f' % elp
			print '\tstart metric calculation'
			for frame in Args['frame']:
				L = []
				elp = timeit.default_timer() - start_time
				print '\t\tstart %i bin frame calculation' %  frame
				if Args['use_synblocks'] == True: 
					syn_path = '%s/%s' % (Args['synblocks_path'],Args['synblocks_files'][sm][s[0]])
					SB = psf.readSynBlocks2(syn_path,resolution,frame+1)
				else: SB = False
				#out_stat = '%s/%s_%s.%iframe.metric.stat' % (Args['out_path'],Args['out_names'][sm],file[:-12],frame)
				M = psf.metricCalc(allCon,resolution,frame=frame,synblocks=SB,metric=Args['metric'])
				M.sort(key=lambda x: Order[x[0]])
				elp = timeit.default_timer() - start_time
				print '\t\tmetric for %i bin frame calculation: %.2f sec' % ( frame, elp)
				if Args['statistics'] != 'only':
					out = '%s/%s_%s.%s.%iframe.metric.bedGraph' % (Args['out_path'],Args['out_names'][sm],file[:-13],Args['metric'],frame)
					f = open(out,'w')
					for m in M: 
						L.append(m[3])
						print >> f, '%s\t%i\t%i\t%f' % (m[0],m[1]*resolution,m[2]*resolution-1,m[3])
					f.close()
				else:
					for m in M: L.append(m[3])
				del M
				elp = timeit.default_timer() - start_time
				print '\t\tmetric for %i bin frame calculating: %.2f sec' % ( frame, elp)
				if Args['statistics'] != False:
					print '\t\tstart calulating metric for randomized contacts'
					randCon = psf.randomizeContacts(allCon)
					elp = timeit.default_timer() - start_time
					print '\t\tcontact randomizing: %.2f sec' % elp
					L = []
					M = psf.metricCalc(randCon,resolution,frame=frame,synblocks=SB,metric=Args['metric'])
					del randCon
					M.sort(key=lambda x: Order[x[0]])
					for m in M: L.append(m[3])
					del M
					print '\t\trandom contact metric calculation: %.2f sec' % elp
					ln = 1.*len(L)
					lb = Args['synblocks_files'][sm][s[0]].split('.')
					h = np.histogram(L,bins=25,range=(0.,.5))
					del L
					plt.plot(h[1][1:],h[0]/ln,color=colorList[c],linestyle='--',label='%s-%s_%s_random' % (lb[0][0],lb[1][0],s[2]))
					elp = timeit.default_timer() - start_time
					print '\tend contact analizying for %i bin frame calculating: %.2f sec' % ( frame, elp)
				else: pass
			del allCon
			elp = timeit.default_timer() - start_time
			print '\t%s end contact analizying %.2f' % (file, elp)
		else: pass
		c += 1
if Args['statistics'] != False:
	plt.legend(fontsize=6)
	plt.savefig('%s/%s.%s.stat.png' % (Args['out_path'],Args['out_names'][sm],Args['metric']),dpi=400)