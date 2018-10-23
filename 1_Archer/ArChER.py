import sys
import timeit
import ArChER_func as Af

Args = {
	'contact_path':'','contact_files':'','genome_path':'','chrom_orders':'',
	'remap_path':'','remap_files':'','out_path':'',
	'resolution':10000, 'confidence':3,'max_difference':33,'model':'balanced'
}

lines = sys.stdin.readlines()
for line in lines:
	parse = line.split('#')
	try:
		key = parse[0].strip().split()[0]
		args = parse[0].strip().split()[2:]
	except IndexError: pass
	else:
		if len(args) == 1:
			try: Args[key] = int(args[0])
			except ValueError: Args[key] = args[0]
			except KeyError: pass
		else: Args[key] = args
		print key, '=', args
start_time = timeit.default_timer()

print '\nStep 0: chromosome indexing...'
l2i = []
fname = [Args['genome_path']+Args['chrom_orders'][i] for i in range(2)]
for i in range(2): l2i.append( Af.iChromIndexing(fname[i]) )
elp = timeit.default_timer() - start_time
print '... chromosome indexing total time:', elp

print '\nStep 1: data reading...'
contactList = []
fname = [Args['contact_path']+Args['contact_files'][i] for i in range(2)]
for i in range(2):
	print '\t' + Args['contact_files'][i]
	contactList.append( Af.iReadPercentelizedContact(fname[i],l2i[i]) )
	elp = timeit.default_timer() - start_time
	print '... locus contact reading end time', elp, len(contactList[i])
elp = timeit.default_timer() - start_time
print '... contact reading total time:', elp

print '\nStep 2: Reading mark points...'
MarkPoints = []
rname = [Args['remap_path']+Args['remap_files'][i] for i in range(2)]
for i in range (2):
	print '\t', rname[i],
	MarkPoints.append( Af.iReadingMarkPoints(rname[i],Args['resolution'],l2i[i],l2i[i-1]) )
	elp = timeit.default_timer() - start_time
	print '\ %i mark point readed for %.2f sec' % (len(MarkPoints[i].keys()), elp)
elp = timeit.default_timer() - start_time
print '... mark point reading total time:', elp

print '\nStep 3: start contact comparing...'
fname = [Args['out_path']+Args['contact_files'][i] for i in range(2)]
for i in range(2):
	out_name = '%s.%im.%iC.%iD.allContacts' % (fname[i],Args['model'],Args['confidence'],Args['max_difference'])
	print '\tstart contact comparing...', out_name
	Dif_Contact = Af.iDifferContact(contactList[i], contactList[1-i], MarkPoints[i], Args['resolution'], Args['model'], Args['confidence'], Args['max_difference'],l2i[i-1])
	elp = timeit.default_timer() - start_time
	print '\tend contact comparing', elp
	
	print '\tDiffer contact writing'
	Af.iPrintDifferContact(Dif_Contact, Args['resolution'], l2i[i], out_name)
	print '\t\tdistant contact writing'
	del Dif_Contact
	elp = timeit.default_timer() - start_time
	print '\tend contact writing', elp

print '...end contact comparing', elp
exit()