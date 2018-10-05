import sys
import timeit
import ArChER_func as af

reload(af)

lines = sys.stdin.readlines()
path = lines[0].split('=')[1].strip()
file = lines[1].split('=')[1].strip().split()
Resolution = int(lines[2].split('=')[1].strip())
path_remap = lines[3].split('=')[1].strip()
file_remap = lines[4].split('=')[1].strip().split()
model = int(lines[5].split('=')[1].strip())
confidence = int(lines[6].split('=')[1].strip())
equal = int(lines[7].split('=')[1].strip())
file_out = lines[8].split('=')[1].strip().split()

file_name = [path+'/'+i for i in file]
Locus_contact_list = []

print 'Step 1: data reading...'
start_time = timeit.default_timer()
for i in range(len(file)):
	print '\t' + file[i]
	File = open(file_name[i]+'.initialContacts','r')
	lines = File.readlines()[1:]
	File.close()
	elp = timeit.default_timer() - start_time
	print '\tend file reading', elp
	Chr_Ind = af.ReadChrIndex(file_name[i]+'.chrInd') 
	Locus_contact_list.append( af.ReadPercentelizedContact(lines) )
	L = len(lines)
	del lines
	elp = timeit.default_timer() - start_time
	print '... %i locus contact readed end time' % L, elp, len(Locus_contact_list[i])
elp = timeit.default_timer() - start_time
print '... contact reading total time:', elp

print 'Step 2: Reading mark points...'
file_name = [path_remap+'/'+i for i in file_remap]
MarkPoints = []

for i in range (0,len(file_name)):
	print '\t', file_remap[i],
	File = open(file_name[i],'r')
	lines = File.readlines()[1:]
	File.close()
	MarkPoints.append( af.ReadingMarkPoints(lines,Resolution) )
	L = len(lines)
	del lines
	elp = timeit.default_timer() - start_time
	print '\t%i mark point for %i loci readed' % (L, len(MarkPoints[i].keys())), elp
elp = timeit.default_timer() - start_time
print '... mark point reading total time:', elp

print 'Step 3: start contact comparing...'

for i in range(2):
	print '\tstart contact comparing...', file_out[i]
	Dif_Contact = af.DifferContactShort(Locus_contact_list[i], Locus_contact_list[1-i], MarkPoints[i], Resolution, model, confidence, equal, 1)
	elp = timeit.default_timer() - start_time
	print '\tend contact comparing', elp
	print '\tDiffer contact writing'
	af.PrintDifferContact(Dif_Contact, Resolution, path+'/%s.%im.%iS.%ieql.allContacts' % (file_out[i],model,confidence,equal))
	print '\t\tdistant contact writing'
	del Dif_Contact
	elp = timeit.default_timer() - start_time
	print '\tend contact writing', elp

print '...end contact comparing', elp
exit()