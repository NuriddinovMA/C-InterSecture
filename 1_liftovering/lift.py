import os
import sys
import timeit
import lift_func as lf

Args = {
	'contact_path':'','contact_files':'','genome_path':'','chrom_orders':'',
	'remap_path':'','remap_files':'','out_path':'',
	'resolution':50000,'agg_frame':[150000],'model':'balanced','inter': False,'dups_filter': 'default'
}

lines = sys.stdin.readlines()
for line in lines:
	parse = line.split('#')
	try:
		key = parse[0].strip().split()[0]
		args = parse[0].strip().split()[2:]
	except IndexError: pass
	else:
		if key == 'agg_frame': Args[key] = [ int(s) for s in args ]
		elif len(args) == 1:
			try: Args[key] = int(args[0])
			except ValueError: Args[key] = args[0]
			except KeyError: pass
		else: Args[key] = args
		#print key, '=', args
start_time = timeit.default_timer()

Args['inter'] == lf.boolean(Args['inter'])

if len(Args['agg_frame']) == 0:
	print 'Using default values of agg_frame: 150000 bp'
	Args['agg_frame'] = [150000,150000]
elif len(Args['agg_frame']) == 1: 
	print 'Given agg_frame value is used for both species'
	Args['agg_frame'].append(Args['agg_frame'][0])
elif len(Args['agg_frame']) > 2:
	print 'too much values of agg_frame'
	exit()
else: pass

if Args['dups_filter'] == 'length' or Args['dups_filter'] == 'coverage' or Args['dups_filter'] == 'deviation': pass
else: 
	print 'Invalid duplicate_filter value!'
	exit()

for key in Args: print key, ' = ', Args[key]

try: os.makedirs(Args['out_path'])
except OSError: pass

print '\nStep 0: chromosome indexing...'
l2i = []
fname = [Args['genome_path']+Args['chrom_orders'][i] for i in range(2)]
for i in range(2): l2i.append( lf.ChromIndexing(fname[i]) )
elp = timeit.default_timer() - start_time
print '... chromosome indexing total time:', elp

print '\nStep 1: data reading...'
contactList = []
fname = [Args['contact_path']+Args['contact_files'][i] for i in range(2)]
for i in range(2):
	print '\t' + Args['contact_files'][i]
	contactList.append( lf.iReadPercentelizedContact(fname[i],l2i[i]) )
	elp = timeit.default_timer() - start_time
	print '... locus contact reading end time', elp, len(contactList[i])
elp = timeit.default_timer() - start_time
print '... contact reading total time:', elp

print '\nStep 2: Reading mark points...'
MarkPoints = []
rname = [Args['remap_path']+Args['remap_files'][i] for i in range(2)]
for i in range (2):
	print '\t', rname[i],
	MarkPoints.append( lf.iReadingMarkPoints(rname[i],Args['resolution'],l2i[i],l2i[i-1],Args['agg_frame'][i]) )
	elp = timeit.default_timer() - start_time
	print '\ %i mark point readed for %.2f sec' % (len(MarkPoints[i].keys()), elp)
elp = timeit.default_timer() - start_time
print '... mark point reading total time:', elp

print '\nStep 3: start contact comparing...'
fname = [Args['out_path']+Args['contact_files'][i] for i in range(2)]
for i in range(2):
	out_name = '%s.%s.liftContacts' % (fname[i],Args['model'])
	print '\tstart contact comparing...', out_name
	Dif_Contact = lf.iDifferContact(contactList[i], contactList[1-i], MarkPoints[i], Args['resolution'], Args['inter'], Args['model'], Args['dups_filter'], l2i[i-1],fname[i])
	elp = timeit.default_timer() - start_time
	print '\tend contact comparing', elp
	
	print '\tDiffer contact writing'
	lf.iPrintDifferContact(Dif_Contact[0], Args['resolution'], l2i[i], out_name, False)
	print '\t\tdistant contact writing'
	lf.iPrintDifferContact(Dif_Contact[1], Args['resolution'], l2i[i], out_name, True)
	print '\t\tdropped dups contact writing'
	del Dif_Contact
	elp = timeit.default_timer() - start_time
	print '\tend contact writing', elp

print '...end contact comparing', elp
exit()
