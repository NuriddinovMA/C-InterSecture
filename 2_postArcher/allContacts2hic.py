import os
import sys
import timeit
import postArcher_func as psAf
reload(psAf)

print 'Step 1: Sites Reading...'
start_time = timeit.default_timer()

Args = {
	'contact_path':'','samples':[],'contact_files':[],
	'frame':[], 'loci': [], 'statistics': False,
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
		elif key == 'use_synblocks' or key == 'statistics' or key == 'use_pre' or key == 'use_loci': Args[key] = psAf.boolean(args[0])
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

for key in Args.keys(): print '\t%s =' % key, Args[key]
command = "java -jar %s pre" % Args['path_to_juicer']
print '\tcommand', command

RH = {
	10000: "1000000,200000,100000,50000,20000,10000",
	20000: "1000000,200000,100000,40000,20000",
	25000: "1000000,200000,100000,50000,25000",
	40000: "1000000,200000,80000,40000",
	50000: "1000000,200000,100000,50000"
}
for sm in range(len(Args['samples'])):
	if Args['contact_files'][sm][0] == 'all' or Args['contact_files'][sm][0][:3] == 'key': files = os.listdir( '%s/%s/' % (Args['contact_path'],Args['samples'][sm]))
	else: files = [ a for a in Args['contact_files'][sm] ]
	if Args['contact_files'][sm][0][:3] == 'key':
		for fi in range(len(files)-1,-1,-1):
			if files[fi].find(Args['contact_files'][sm][0][4:]) == -1: del files[fi]
	for file in files:
		s = file.split('.')
		fname = '%s/%s/%s' % (Args['contact_path'],Args['samples'][sm],file)
		out = '%s/%s_%s' % (Args['out_path'],Args['out_names'][sm],file[:-12])
		print 'read file', fname
		print 'out file', out
		if Args['chrom_sizes'].has_key(s[0]) == True and s[-1] == 'allContacts':
			resolution = int(s[2][:-2])*1000
			G = '%s/%s' % (Args['chrom_path'],Args['chrom_sizes'][s[0]])
			if Args['use_pre'] == True:
				try:
					f = open(out + 'Reference.pre')
					f.close()
					f = open(out + 'Query.pre')
					f.close()
				except IOError:
					print '\tpre not found!'
					Args['use_pre'] = False
			if Args['use_pre'] == False:
				print '\tstart reading', s[0], s[1], s[2], elp
				Order = psAf.ChromIndexing(G)
				allCon = psAf.readContacts(fname,Order,resolution)
				print '\tgenerate pre', out
				M = psAf.JuiceboxPre(allCon,Order,resolution,out)
				del allCon
				print '\tpre writing', elp
			R = RH[resolution]
			F = out + '.Reference.pre'
			O = out + '.Reference.hic'
			os.system( command + " " + F + " " + O + " " + G + " " + "-r" + " " + R + " "+ "-n")
			F = out + '.Query.pre'
			O = out + '.Query.hic'
			os.system( command + " " + F + " " + O + " " + G + " " + "-r" + " " + R + " "+ "-n")
		else: pass