import pandas as pd

def read1(path):
	A = []
	start_time = timeit.default_timer()
	f = open(path,'r')
	lines = f.readlines()
	f.close()
	print timeit.default_timer() - start_time
	for i in range(len(lines)-1,-1,-1):
		parse = lines[i].split()
		del lines[i]
		A.append( ( int(parse[1][3:]),int(parse[2])/10000,int(parse[4][3:]),int(parse[5]/10000) ) )
		if ((i+1) % 1000000) == 0: print timeit.default_timer() - start_time
	return A

def read2(path):
	A = []
	start_time = timeit.default_timer()
	f = open(path,'r')
	i = 0
	while True:
		i += 1
		parse = f.readline().split()
		A.append( ( int(parse[1][3:]),int(parse[2])/10000,int(parse[4][3:]),int(parse[5]/10000) ) )
		if i % 1000000 == 0: print timeit.default_timer() - start_time
	f.close()
	return A

def read_pandas(path):
	A = pd.read_csv(path, sep='\t')
	print timeit.default_timer() - start_time
	return A
