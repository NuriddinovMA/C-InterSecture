import os
import sys
import timeit
import lift_func as lf

Args = {
	'track_path':'','track_files':[],'genome_path':'','chrom_orders':[],
	'remap_path':'','remap_files':[],'out_path':''
}

lines = sys.stdin.readlines()
for line in lines:
	parse = line.split('#')
	try:
		key = parse[0].strip().split()[0]
		args = parse[0].strip().split()[2:]
	except IndexError: pass
	else:
		try:
			if key == 'chrom_orders': Args[key].append(args)
			elif key == 'track_files' or key == 'remap_files': Args[key] = args
			else: Args[key] = args[0]
		except KeyError: pass
start_time = timeit.default_timer()

for key in Args: print key, ' = ', Args[key]

try: os.makedirs(Args['out_path'])
except OSError: pass

num = len(Args['track_files'])

for i in range(num):
	print '\nStep %i.0: chromosome indexing...' % i
	l2i = []
	fname = [Args['genome_path']+Args['chrom_orders'][i][j] for j in range(2)]
	for j in range(2): l2i.append( lf.ChromIndexing(fname[j]) )
	elp = timeit.default_timer() - start_time
	print '... chromosome indexing total time:', elp

	print '\nStep %i.1: Reading mark points...' % i
	rname = Args['remap_path']+Args['remap_files'][i]
	MarkPoints = lf.lftReadingMarkPoints(rname,l2i[0],l2i[1])
	elp = timeit.default_timer() - start_time
	print '... mark point reading total time:', elp

	print '\nStep %i.2: start contact lifting over...' % i
	fname = Args['track_path']+Args['track_files'][i]
	out_name = Args['out_path']+Args['track_files'][i][:-3] + Args['remap_files'][i][:-4] + 'bed'
	lf.lftRough(MarkPoints, l2i[0],l2i[1], fname, out_name)
	elp = timeit.default_timer() - start_time
	print '\tend contact writing', elp

print '...end contact lifting over', elp
exit()