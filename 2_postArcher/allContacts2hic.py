import os
import timeit
import postArcher_func as pAf
reload(pAf)

print 'Step 1: Sites Reading...'
start_time = timeit.default_timer()

resolution = 40000
base_dir = '//fishminja/share/ArChER/'
file_path = 'HumanMouse/','HumanChicken/','MouseChicken/'
start_time = timeit.default_timer()
order_path = 'chrSizes/'
order_file = {'IMR90':'hg38','MKEF':'mm10','CEF':'galGal5'}
path_to_juicer = "C:/Desktop/juicebox/juicertools.jar"
command = "java -jar " + path_to_juicer + ' pre'
frame = 8
for p in range(len(file_path)):
	files = os.listdir(base_dir+file_path[p])
	elp = timeit.default_timer() - start_time
	print file_path[p][:-1], elp
	for file in files:
		s = file.split('.')
		fname = base_dir+file_path[p]+file
		if s[-1] == 'allContacts':
			resolution = int(s[4][:-2])*1000
			G = base_dir+order_path+order_file[s[0]] +'.chr.sizes'
			try:
				f = open(fname+'.Reference.pre')
				f.close()
				f = open(fname+'.Query.pre')
				f.close()
			except IOError:
				print '\tstart reading', s[0], s[4], elp
				Order = pAf.chrSortOrder(G)
				allCon = pAf.readContacts2(fname,Order,resolution)
				M = pAf.JuiceboxPre(allCon,Order,resolution,fname)
				del allCon
				print '\tpre writing', elp
			R = "1000000,200000,40000"
			if resolution < 40000: R += "," + str(resolution)
			F = fname + '.Reference.pre'
			O = fname + '.Reference.hic'
			os.system( command + " " + F + " " + O + " " + G + " " + "-r" + " " + R + " "+ "-n")
			F = fname + '.Query.pre'
			O = fname + '.Query.hic'
			os.system( command + " " + F + " " + O + " " + G + " " + "-r" + " " + R + " "+ "-n")
		else: pass