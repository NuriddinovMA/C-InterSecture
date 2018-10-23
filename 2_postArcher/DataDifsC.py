import timeit
import ArcherStatFunction as FS
reload(FS)

print 'Step 3: Sites Reading...'
start_time = timeit.default_timer()

Resolution = 40000
base_dir = 'C:/Desktop/Arbeit/AllData/ArChER/differ/'
file_path = 'HM/'
suf = '.33N.1C.40kb.5Mb.0m.3S.33eql'
sp_name = ['IMR90',]#'MKEF')

hS = []
difCon = []
eqlCon = []
allCon = []

start_time = timeit.default_timer()


for i in range (len(sp_name)):
	con_name = base_dir+file_path+sp_name[i]+suf 
	elp = timeit.default_timer() - start_time
	print '\tstart %s unmapped contact stat total time:' % sp_name[i], elp
	File = open(con_name+'.allContacts','r')
	lines = File.readlines()[1:]
	File.close()
	hS = FS.ReadingBins(lines)
	del lines
	allCon.append( FS.BinsStat(hS, Resolution, 'c', 0) )
	hS = []
	elp = timeit.default_timer() - start_time
	print '\t...differ contact stat total time:', elp
'''for i in range (2):
	File = open(dif_name[i]+suf2+'.equalContacts','r')
	lines = File.readlines()[1:]
	File.close()
	hS = FS.ReadingBins(lines)
	del lines
	eqlCon.append( FS.BinsStat(hS, Resolution, 'c', 0) )
	hS = []
	elp = timeit.default_timer() - start_time
	print '\t...equal contacts readed', dif_name[i], elp'''

