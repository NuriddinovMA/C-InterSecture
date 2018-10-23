import os
import timeit
import postArcher_func as pAf
reload(pAf)

print 'Step 1: Sites Reading...'
start_time = timeit.default_timer()

resolution = 40000
base_dir = '//fishminja/share/ArChER/'
file_path = 'HumanMouse/','HumanChicken/','MouseChicken/'
out_path = 'tracks/'
start_time = timeit.default_timer()
order_path = 'chrSizes/'
order_file = {'IMR90':'hg38','MKEF':'mm10','CEF':'galGal5'}
frame = 8
for p in range(1):#len(file_path)):
	files = os.listdir(base_dir+file_path[p])
	elp = timeit.default_timer() - start_time
	print file_path[p][:-1], elp
	for file in files[:-2]:
		s = file.split('.')
		if s[-1] == 'allContacts' and s[1] == 'prc':
			print '\tstart reading', s[0], s[4], elp
			resolution = int(s[4][:-2])*1000
			Order = pAf.chrSortOrder(base_dir+order_path+order_file[s[0]] +'.chr.sizes')
			allCon = pAf.readContacts2(base_dir+file_path[p]+file,Order,resolution)
			M = pAf.metricCalc(allCon,resolution,frame=frame)
			del allCon
			out = file_path[p][:-1] + '_%s.%s.%s.%iframe.metric2.bedGraph' % (s[0],s[1],s[4],frame)
			M.sort(key=lambda x: Order[x[0]])
			elp = timeit.default_timer() - start_time
			print '\tmetric calculation', elp
			f = open(base_dir+out_path+out,'w')
			for m in M: print >> f, '%s\t%i\t%i\t%f' % (m[0],m[1]*resolution,m[2]*resolution-1,m[3])
			f.close()
			del M
			print '\tmetric writing', elp
		else: pass