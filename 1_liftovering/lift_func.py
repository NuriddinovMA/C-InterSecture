import math
import numpy as np
import timeit
import sys

def colorList():
	color = []
	for i in range(0,256,85):
		for j in range(0,256,85):
			for k in range(0,256,85): color.append('%i,%i,%i' % (i,j,k))
	return color

def HashTry(Hash, key):
	t = 1
	try:
		Hash[key]
	except KeyError:
		t = 0
	return t

def FindKey(Hash, key):
	if key != 'interchromosome':
		for k in range(key,-1,-1):
			try: 
				Hash[k]
				return k
			except KeyError: pass
 	else: return key

def boolean(x):
	if x == 'False' or x == 'false' or x == 'f' or x == 'F': return False
	elif x == 'True' or x == 'true' or x == 't' or x == 'T': return True
	else: return x

def ChromIndexing(path):
	ChrInd = {}
	f = open(path,'r')
	lines = f.readlines()
	f.close()
	for i in range(len(lines)):
		parse = lines[i].split()
		try: 
			ChrInd[parse[0]] = i+1
			ChrInd[i+1] = parse[0]
		except IndexError: break
	del lines
	return ChrInd

def iReadInitialContact(path,ChrIdxs): #Reading Contact from database file
	contactHash = {}
	start_time = timeit.default_timer()
	f = open(path+'.initialContacts','r')
	lines = f.readlines()
	f.close()
	elp = timeit.default_timer() - start_time
	ln = len(lines)
	print '\t\t readed contact count: %i, time elapsed: %.2f' % (ln,elp)
	for i in range(ln-1,0,-1):
		parse = lines[i].split('\t')
		del lines[i]
		try:
			contactHash[ChrIdxs[parse[0]], int(parse[1]), ChrIdxs[parse[2]], int(parse[3])] = int(parse[4]),int(parse[5]),int(parse[6]),int(parse[7]),int(parse[8]),float(parse[9])
			if (ln-i) % 10000000 == 0:
				elp = timeit.default_timer() - start_time
				print '\t\t contact reading: %i, time elapsed: %.2f, memory sized: %.2fMb' % (ln-i,elp,sys.getsizeof(contactHash)/1024./1024.)
		except KeyError: pass
	return contactHash

def readPSHash(path,resolution):
	psHash = {}
	f = open(path,'r')
	lines = f.readlines()
	f.close()
	for i in range(len(lines)):
		parse = np.float32(lines[i].split()).tolist()
		if len(parse[2:]) == 0: 
			psHash = False
			break
		else: psHash[parse[0]] = parse[2:]
	return psHash

def iPercentelizedContactStatistic(path): #Reading Contact from database file
	Stat = [[0 for j in range(10)] for i in range(5)], [[0 for j in range(10)] for i in range(5)], [[0 for j in range(10)] for i in range(5)]
	start_time = timeit.default_timer()
	f = open(path,'r')
	lines = f.readlines()
	f.close()
	ln = len(lines)
	for i in range(ln-1,-1,-1):
		parse = lines[i].split()
		del lines[i]
		if parse[0] == parse[2]:
			k1 = int(parse[1])
			k2 = int(parse[3])
			c = int(parse[4])
			d1 = int(parse[5])
			d2 = int(parse[6])
			d = d2-d1
			f = math.fabs( ( d1+d2 ) / 2.0 - c) 
			if math.fabs(k2-k1) == 2:
				Stat[0][0][(c-1)/10] += 1
				Stat[1][0][(c-1)/10] += d1
			elif math.fabs(k2-k1) == 5:
				Stat[0][1][(c-1)/10] += 1
				Stat[1][1][(c-1)/10] += d1
			elif math.fabs(k2-k1) == 10:
				Stat[0][2][(c-1)/10] += 1
				Stat[1][2][(c-1)/10] += d1
			elif math.fabs(k2-k1) == 20:
				Stat[0][3][(c-1)/10] += 1
				Stat[1][3][(c-1)/10] += d1
			elif math.fabs(k2-k1) == 50:
				Stat[0][4][(c-1)/10] += 1
				Stat[1][4][(c-1)/10] += d1
			else: pass
		if (ln-i) % 10000000 == 0:
			elp = timeit.default_timer() - start_time
			print '\t\tLocus contact reading:', ln-i, 'time elapsed:', elp, 'memory sized:', sys.getsizeof(Stat)
	for i in range(5):
		for j in range(10):
			Stat[2][i][j] = Stat[1][i][j]/Stat[0][i][j]
	return Stat

def netParser(name,**kwargs):
	try: min_length = kwargs['min_length']
	except KeyError: min_length = 0
	try: gap_length = kwargs['gap_length']
	except KeyError: gap_length = 100000
	f = open(name, 'r')
	lines = f.readlines()
	f.close()
	parsedNet = {},{}
	id = 0
	for line in lines:
		if line[0] == 'n': 
			key = line.split()[1]
			parsedNet[0][key] = []
		else:
			parse = line.split()
			start1,ln1,start2,ln2 = int(parse[1]),int(parse[2]),int(parse[5]),int(parse[6])
			dir = int(parse[4]+'1')
			fi1,fi2 = start1+ln1,start2+ln2
			if parse[0] == 'fill':
				id += 1
				data = key,start1,fi1,parse[3],start2,fi2,dir,id
				parsedNet[0][key].append(data)
				parsedNet[1][id] = []
			else:
				data = key,start1,fi1,parse[3],start2,fi2,dir,id
				parsedNet[1][id].append(data)
	del lines

	for key in parsedNet[0]:
		parsedNet[0][key].sort()
		pN = parsedNet[0][key]
		
		for i in range(len(pN)):
			pNj = parsedNet[1][pN[i][-1]]
			chrm1,st1,fi1,chrm2,st2,fi2,dir,id = pN[i][0],pN[i][1],pN[i][2],pN[i][3],pN[i][4],pN[i][5],pN[i][6],pN[i][7]
			if dir > 0:
				for j in range(len(pNj)):
					st1h,st2h = pNj[j][2],pNj[j][5]
					pNj[j] = [chrm1,st1,pNj[j][1],chrm2,st2,pNj[j][4],dir,id]
					st1,st2 = st1h,st2h
					# pNj.append((chrm1,st1,fi1,chrm2,st2,fi2,dir,id))
				pNj.append([chrm1,st1,fi1,chrm2,st2,fi2,dir,id])
			else:
				for j in range(len(pNj)):
					st1h,fi2h = pNj[j][2],pNj[j][4]
					pNj[j] = [chrm1,st1,pNj[j][1],chrm2,pNj[j][5],fi2,dir,id]
					st1,fi2 = st1h,fi2h
					# pNj.append((chrm1,st1,fi1,chrm2,st2,fi2,dir,id))
				pNj.append([chrm1,st1,fi1,chrm2,st2,fi2,dir,id])
		
		for i in range(len(pN)-1):
			try: parsedNet[1][pN[i][7]]
			except KeyError: pass
			else:
				for j in range(i+1,len(pN)):
					try: parsedNet[1][pN[j][7]]
					except KeyError: pass
					else:
						gap0 = pN[j][1] - pN[i][2]
						if pN[i][6] > 0: gap1 = pN[j][4] - pN[i][5]
						else: gap1 = pN[i][4] - pN[j][5]
						
						if pN[i][3] != pN[j][3] or pN[i][6] != pN[j][6]: pass
						elif pN[i][6] > 0 and pN[i][5] > pN[j][4]: pass
						elif pN[i][6] < 0 and pN[i][4] < pN[j][5]: pass
						elif abs(gap0 - gap1) > gap_length or abs(math.log(1.*(gap0+1)/(gap1+1),2))>1.5: pass
						else:
							# parsedNet[0][key][i] = (pN[i][0],pN[j][1],pN[i][2],pN[i][3],min(pN[i][4],pN[j][4]),max(pN[i][5],pN[j][5]),pN[i][6],pN[i][7])
							# pN[i] = (pN[i][0],pN[j][1],pN[i][2],pN[i][3],min(pN[i][4],pN[j][4]),max(pN[i][5],pN[j][5]),pN[i][6],pN[i][7])
							parsedNet[1][pN[i][7]].extend(parsedNet[1][pN[j][7]])
							del parsedNet[1][pN[j][7]]

		for i in range(len(pN)-1,-1,-1):
			try: parsedNet[1][pN[i][7]]
			except KeyError: del parsedNet[0][key][i]
			else:
				if parsedNet[1][pN[i][7]][0][2] - parsedNet[1][pN[i][7]][0][1] < min_length:
					del parsedNet[1][pN[i][7]]
					del parsedNet[0][key][i]
				else:
					for j in parsedNet[1][pN[i][7]]: j[-1] = pN[i][7]
	pN = []
	Keys = parsedNet[1].keys()
	for key in Keys: 
		pN.append( parsedNet[1][key] )
		del parsedNet[1][key]
	del Keys
	del parsedNet
	pN.sort()
	return pN

def net2pre(parsedNet,out):
	markPoints = []
	f1 = open(out+'.pre.mark', 'w')
	f2 = open(out+'.2D.ann', 'w')
	
	print >> f1, 'chr1\tstart1\tend1\tchr2\tstart2\tend2\tid'
	print >> f2, 'chr1\tstart1\tend1\tchr2\tstart2\tend2\tcolor\tcomment'
	for i in parsedNet:
		a = 0
		if i[0][6] > 1: st,fi = i[0][4],i[-1][5]
		else: st,fi = i[-1][4],i[0][5]
		print >> f2, '%s\t%i\t%i\t%s\t%i\t%i\t%s\t%s:%i-%i:%i' % (i[0][0],i[0][1],i[-1][2],i[0][0],i[0][1],i[-1][2],'0,255,0',i[0][3],st,fi,i[0][6])
		for j in i:
			try: print >> f1, '%s\t%i\t%i\t%s\t%i\t%i\t%i\t%i' % (j[0],j[1],j[2],j[3],j[4],j[5],j[6],j[7])
			except TypeError: a = 1
		if a == 1: print i
	f1.close()
	f2.close()
	
	f3 = open(out+'.mark', 'w')
	for i in parsedNet:
		dir = i[0][6]
		if dir > 0:
			for j in range(len(i)):
				ln1,ln2 = i[j][2]-i[j][1],i[j][5]-i[j][4]
				if ln1 > 300: 
					k = ln1/150.0
					sp1,sp2 = ln1/k,ln2/k
					h = []
					for n in range(int(k)): h.append((i[j][1]+sp1*n,i[j][4]+sp2*n))
					h.append( (i[j][2],i[j][5] ) )
					for n in range(int(k)): markPoints.append( (i[j][0],int(h[n][0]),int(h[n+1][0]),i[j][3],int(h[n][1]),int(h[n+1][1]),i[j][6],i[j][7]) )
				else: markPoints.append(i[j])
				try:
					gap1,gap2 = i[j+1][1]-i[j][2],i[j+1][4]-i[j][5]
					if gap1 > 300:
						k = gap1/150.0
						sp1,sp2 = gap1/k,gap2/k
						h = []
						for n in range(int(k)): h.append((i[j][2]+sp1*n,i[j][5]+sp2*n))
						h.append( (i[j+1][1],i[j+1][4] ) )
						for n in range(int(k)): markPoints.append( (i[j][0],int(h[n][0]),int(h[n+1][0]),i[j][3],int(h[n][1]),int(h[n+1][1]),i[j][6],str(i[j][7])+'_gap') )
					else: markPoints.append( (i[j][0],i[j][2],i[j+1][1],i[j][3],i[j][5],i[j+1][4],i[j][6],str(i[j][7])+'_gap') )
				except IndexError: pass
		else:
			for j in range(len(i)):
				ln1,ln2 = i[j][2]-i[j][1],i[j][5]-i[j][4]
				if ln1 > 300: 
					k = ln1/150.0
					sp1,sp2 = ln1/k,ln2/k
					h = []
					for n in range(int(k)): h.append((i[j][1]+sp1*n,i[j][5]-sp2*n))
					h.append( (i[j][2],i[j][4] ) )
					for n in range(int(k)): markPoints.append( (i[j][0],int(h[n][0]),int(h[n+1][0]),i[j][3],int(h[n+1][1]),int(h[n][1]),i[j][6],i[j][7]) )
				else: markPoints.append(i[j])
				try:
					gap1,gap2 = i[j+1][1]-i[j][2],i[j][4]-i[j+1][5]
					if ln1 > 300: 
						k = ln1/150.0
						sp1,sp2 = ln1/k,ln2/k
						h = []
						for n in range(int(k)): h.append((i[j][2]+sp1*n,i[j][4]-sp2*n))
						h.append( (i[j+1][1],i[j+1][5] ) )
						for n in range(int(k)): markPoints.append( (i[j][0],int(h[n][0]),int(h[n+1][0]),i[j][3],int(h[n+1][1]),int(h[n][1]),i[j][6],str(i[j][7])+'_gap') )
					else: markPoints.append( (i[j][0],i[j][2],i[j+1][1],i[j][3],i[j+1][5],i[j][4],i[j][6],str(i[j][7])+'_gap') )
				except IndexError: pass

	for i in markPoints: print >> f3,'%s\t%i\t%i\t%s\t%i\t%i\t%i\t%s' % (i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7])
	f3.close()

def lftReadingMarkPoints(path, ChrIdxs1, ChrIdxs2):  #Creating Mark Point List to convert MarkPoint from species to species 
	ObjCoorMP = {}
	f = open(path, 'r')
	lines = f.readlines()
	f.close()
	for i in range(len(lines)-1,0,-1):
		try:
			parse = lines[i].split()
			del lines[i]
			try: N1 = ChrIdxs1[parse[0]]
			except KeyError: N1 = ChrIdxs1[parse[0][3:]]
			c1 = ( int(parse[1]) + int(parse[2]) ) / 400
			try: N2 = ChrIdxs2[parse[3]]
			except KeyError: N1 = ChrIdxs1[parse[3][3:]]
			c1 = ( int(parse[1]) + int(parse[2]) ) / 400
			c2 = ( int(parse[4]) + int(parse[5]) ) / 400
			if HashTry(ObjCoorMP, (N1,c1) ) == 0: ObjCoorMP[N1,c1] = set([ (N2,c2) ])
			else: ObjCoorMP[N1,c1].add( (N2,c2) )
		except IndexError: print 'Exception: IndexError',line
		except KeyError: pass
	return ObjCoorMP

def lftRough(ObjCoorMP, ChrIdxs1, ChrIdxs2, bed, out):
	lift = []
	f = open(bed, 'r')
	lines = f.readlines()
	f.close()
	f = open(out, 'w')
	for i in range(len(lines)):
		parse = lines[i].split()
		try: N1 = ChrIdxs1[parse[0]]
		except KeyError: N1 = ChrIdxs1[parse[0][3:]]
		except IndexError: parse
		c1 = int(parse[1]) / 200
		c2 = int(parse[2]) / 200
		for j in range(c1,c2+1,1):
			try:
				lifted = ObjCoorMP[N1,j]
				for l in lifted: print >> f, ChrIdxs2[l[0]], l[1]*200, l[1]*200+199, parse[3]
			except KeyError: pass
	f.close()
	return lift

def iReadingMarkPoints(path, resolution, ChrIdxs1,ChrIdxs2,agg):  #Creating Mark Point List to convert MarkPoint from species to species 
	ObjCoorMP = {},{},{}
	key = ''
	frame = agg/resolution
	f = open(path+'.mark', 'r')
	lines = f.readlines()
	f.close()
	for i in range(len(lines)-1,0,-1):
		try:
			parse = lines[i].split()
			del lines[i]
			try: N1 = ChrIdxs1[parse[0]]
			except KeyError: N1 = ChrIdxs1[parse[0][3:]]
			c1 = ( int(parse[1]) + int(parse[2]) ) / 100 * 50
			try: N2 = ChrIdxs2[parse[3]]
			except KeyError: N2 = ChrIdxs2[parse[3][3:]]
			b2 = ( int(parse[4]) + int(parse[5]) ) / (2*resolution)
			if HashTry(ObjCoorMP[1], (N2,b2) ) == 0: ObjCoorMP[1][N2,b2] = set([ (N1,c1) ])
			else: ObjCoorMP[1][N2,b2].add( (N1,c1) )
			if HashTry(ObjCoorMP[2], (N1,c1) ) == 0: ObjCoorMP[2][N1,c1] = set([ (N2,b2) ])
			else: ObjCoorMP[2][N1,c1].add( (N2,b2) )
		except IndexError:  print 'Exception: IndexError',line
		except KeyError: pass
	for i in ObjCoorMP[1]:
		for j in ObjCoorMP[1][i]:
			k = 1.0/len(ObjCoorMP[1][i])
			key = j[0], j[1]/resolution
			if HashTry(ObjCoorMP[0], key ) == 0:
				ObjCoorMP[0][key] = [{}]
				ObjCoorMP[0][key][0][i] = [1,k]
			else:
				p = 0
				for n in range(len(ObjCoorMP[0][key])):
					if (HashTry(ObjCoorMP[0][key][n], i) == 1):
						ObjCoorMP[0][key][n][i] = [ObjCoorMP[0][key][n][i][0] + 1, ObjCoorMP[0][key][n][i][1] + k]
						sorted(ObjCoorMP[0][key][n][i])
						p = 1
						break
					elif ( (ObjCoorMP[0][key][n].keys()[0][0] == i[0]) and (math.fabs(ObjCoorMP[0][key][n].keys()[0][1] - i[1]) <= frame) ):
						ObjCoorMP[0][key][n][i] = [1, k]
						sorted(ObjCoorMP[0][key][n][i])
						p = 1
						break
					else: pass
				if p == 0: ObjCoorMP[0][key].append({i:[1,k]})
	N = [0,0,0]
	headKeys = ObjCoorMP[0].keys()
	for key in headKeys:
		S = 0.0
		for n in ObjCoorMP[0][key]: 
			for k in n: S += n[k][0]
		for n in range(len(ObjCoorMP[0][key])-1,-1,-1):
			lowKeys = ObjCoorMP[0][key][n].keys()
			for k in lowKeys:
				ObjCoorMP[0][key][n][k] = ( round( ObjCoorMP[0][key][n][k][0]/S, 2), round( ObjCoorMP[0][key][n][k][1], 2) )
				if (ObjCoorMP[0][key][n][k][0] == 0.0) or (ObjCoorMP[0][key][n][k][1] == 0.0):
					N[0] += 1
					del ObjCoorMP[0][key][n][k]
			if len(ObjCoorMP[0][key][n]) == 0:
				N[1] += 1
				del ObjCoorMP[0][key][n]
		if len(ObjCoorMP[0][key]) == 0:
			N[2] += 1
			del ObjCoorMP[0][key]
	print '\t\tdropped remapped bins', N
	sorted(ObjCoorMP[0])
	#return ObjCoorMP
	return ObjCoorMP[0]

def iDuplicateContact(Contact_disp_L, ObjCoorMP1, ObjCoorMP2):
	c = [0,0,0,0,0,0,0]
	c[-1] = ObjCoorMP1[0]*ObjCoorMP2[0]*ObjCoorMP1[1]*ObjCoorMP2[1]
	c[0] = c[-1]*Contact_disp_L[0]
	c[1] = c[-1]*Contact_disp_L[1]
	c[2] = c[-1]*Contact_disp_L[2]
	c[3] = c[-1]*Contact_disp_L[3]
	c[4] = c[-1]*Contact_disp_L[4]
	c[-2] = c[-1]*Contact_disp_L[-1]
	return c
	
def iLabel(KeyList, resolution, ChrIdxs2):
	lb = ''
	if len(KeyList) < 3:
		KeyList = list(KeyList)
		KeyList.sort()
		for k in KeyList:  
			lb += '%s:%i-' % (ChrIdxs2[k[0]],k[1]*resolution)
	else:
		KeyList = list(KeyList)
		KeyList.sort()
		lb = '%s:%i-' % (ChrIdxs2[KeyList[0][0]],KeyList[0][1]*resolution)
	if lb == '': lb = '-------'
	return lb

def iChouseBest(writed, dupled, crit):
	result = True
	if crit == 'coverage':
		if writed[-2]*writed[-3] < dupled[-2]*dupled[-3] : results = True
		else: results = False
	elif crit == 'deviation':
		if writed[5] < dupled[5] : results = True
		else: results = False
	elif crit == 'length':
		if (writed[-2] < dupled[-2]) or (dupled[-2] < 0): results = True
		else: results = False
	else:
		if (writed[-1] < dupled[-1]): results = True
		else: results = False
	
	return result

def iDifferContactStat(Contact_disp_0, Contact_disp_1, ObjCoorMP, model):
	DifferContact = {},{}
	Stat = [ [0 for i in range(101)] for j in range(3)]
	for i in Contact_disp_0:
		key1 = i[:2]
		key2 = i[2:]
		Stat[0][int(Contact_disp_0[i][0])] += 1
		if  (HashTry(ObjCoorMP, key1) == 1) and (HashTry(ObjCoorMP, key2) == 1):
			end1 = len(ObjCoorMP[key1])
			end2 = len(ObjCoorMP[key2])
			k = 0
			for j1 in range( end1):
				for j2 in range( end2):
					for k1 in ObjCoorMP[key1][j1]:
						for k2 in ObjCoorMP[key2][j2]:
							if HashTry(Contact_disp_1,k1+k2) == 1 or HashTry(Contact_disp_1,k2+k1) == 1: k += 1
							else: pass
			if k == 0: Stat[2][int(Contact_disp_0[i][0])] += 1
			else: pass
		else: Stat[1][int(Contact_disp_0[i][0])] += 1
	for i in range(101):
		Stat[1][i] = Stat[0][i] - Stat[1][i]
		Stat[2][i] = Stat[1][i] - Stat[2][i]
	return Stat

def iDifferContact(Contact_disp_0, Contact_disp_1, ObjCoorMP, resolution, inter, model, criteria, ChrIdxs2,stat_out):
	DifferContact = {}
	Dups = {}
	Statistic = ['all','remappable','processed','unique','duplicated','dropped'],[0,0,set([]),set([]),[],0],[set([]),set([]),set([]),set([]),set([]),set([])]
	for i in Contact_disp_0:
		key1 = i[:2]
		key2 = i[2:]
		Statistic[1][0] += 1
		Statistic[2][0] |= set([key1,key2,])
		if (key1[0] != key2[0]) and inter == False: pass
		elif (HashTry(ObjCoorMP, key1) == 1) and (HashTry(ObjCoorMP, key2) == 1):
			end1 = len(ObjCoorMP[key1])
			end2 = len(ObjCoorMP[key2])
			k = 0
			c = [0,0,0,0,0,0,0,set(),set()]
			Statistic[1][1] += 1
			Statistic[2][1] |= set([key1,key2,])
			if end1 == 0: print "remapping error!!!", ObjCoorMP[key1]
			elif end1 == 1: Statistic[2][3] |= set([key1,])
			else: Statistic[2][4] |= set([key1,])
			if end2 == 0: print "remapping error!!!", ObjCoorMP[key2]
			elif end2 == 1: Statistic[2][3] |= set([key2,])
			else: Statistic[2][4] |= set([key2,])
			for j1 in range(end1):
				for j2 in range(end2):
					c = [0,0,0,0,0,0,0,set(),set()]
					k = 0
					for k1 in ObjCoorMP[key1][j1]:
						for k2 in ObjCoorMP[key2][j2]:
							if (k1[0] != k2[0]) and (inter == False): pass
							elif HashTry(Contact_disp_1,k1+k2) == 1:
								dc = iDuplicateContact(Contact_disp_1[k1+k2], ObjCoorMP[key1][j1][k1], ObjCoorMP[key2][j2][k2])
								c[0] += dc[0]
								c[1] += dc[1]
								c[2] += dc[2]
								c[3] += dc[3]
								c[4] += dc[4]
								c[-4] += dc[-2]
								c[-3] += dc[-1]
								c[-2].add(k1)
								c[-1].add(k2)
								k += 1
							elif HashTry(Contact_disp_1,k2+k1) == 1:
								dc = iDuplicateContact(Contact_disp_1[k2+k1], ObjCoorMP[key1][j1][k1], ObjCoorMP[key2][j2][k2])
								c[0] += dc[0]
								c[1] += dc[1]
								c[2] += dc[2]
								c[3] += dc[3]
								c[4] += dc[4]
								c[-4] += dc[-2]
								c[-3] += dc[-1]
								c[-2].add(k2)
								c[-1].add(k1)
								k += 1
							else: Statistic[2][5] |= set([key1,key2,])
					if k == 0: pass
					else:
						#Statistic[1][2] += 1
						Statistic[1][2] |= set([i,])
						Statistic[2][2] |= set([key1,key2,])
						if model != 'balanced': norm = 1
						else: norm = c[-3]
						if c[-3] != 0:
							cc = (
								int(round((c[0])/norm)), 
								int(round((c[1])/norm)), 
								int(round((c[2])/norm)),
								int(round((c[3])/norm)),
								int(round((c[4])/norm)),
								float(round((c[-4])/norm,2)),
								c[-3],
								set(c[-2]),
								set(c[-1])
							)
							disp1 = max(( Contact_disp_0[i][2]-Contact_disp_0[i][0]),(Contact_disp_0[i][0]-Contact_disp_0[i][1])) 
							disp2 = max((cc[2]-cc[0]),(cc[0]-cc[1]))
							to_write = (iLabel(cc[-2],resolution,ChrIdxs2)[:-1], iLabel(cc[-1],resolution,ChrIdxs2)[:-1], Contact_disp_0[i][0], cc[0], disp1, disp2, Contact_disp_0[i][3], Contact_disp_0[i][4], cc[3], cc[4], cc[-4],cc[-3])
							if HashTry(DifferContact,i) == 0: 
								DifferContact[i] = to_write
								Statistic[1][3] |= set([i])
							else:
								Statistic[1][4].append(i)
								#Statistic[2][4] += [key1,key2,] 
								#if DifferContact[i][-1] < cc[-3] or cc[-3] <= 0 : pass
								if iChouseBest(DifferContact[i],to_write,criteria) == True: 
									if HashTry(Dups,i) == 0: Dups[i] = [to_write,]
									else: Dups[i].append(to_write)
								else:
									if HashTry(Dups,i) == 0: Dups[i] = [DifferContact[i],]
									else: Dups[i].append(DifferContact[i])
									DifferContact[i] = to_write
						else: pass
		else: pass
	Statistic[1][2] = len(Statistic[1][2])
	Statistic[1][5] = len(Statistic[1][4])
	Statistic[1][4] = len( set(Statistic[1][4]) )
	Statistic[1][3] = len( Statistic[1][3] ) - Statistic[1][4]
	Statistic[2][5] = Statistic[2][5] - Statistic[2][2]
	Statistic[2][3] = Statistic[2][3] - Statistic[2][5]
	Statistic[2][4] = Statistic[2][4] - Statistic[2][5]
	f = open(stat_out+'.stat','w')
	for i in range(6): 
		Statistic[2][i] = len(Statistic[2][i])
		print >> f, Statistic[0][i],Statistic[1][i],Statistic[2][i]
	del Statistic
	return DifferContact,Dups

def iPrintDifferContact(data, resolution, ChrIdxs, out, dups,**kwargs):
	if dups: out += '.dups'
	try: psHash=kwargs['scoring']
	except KeyError: psHash == False
	f = open(out, 'w')
	Keys = sorted(data.keys(),key=lambda x: (x[0],x[2],x[1],x[3]))
	print >> f, 'chr1_observed\tpos1_observed\tchr2_observed\tpos2_observed\tremap1_control\tremap2_control\tobserved_contacts\tcontrol_contacts\tobserved_deviations\tcontrol_deviations\tobserved_coverages_pos1\tobserved_coverages_pos2\tcontrol_coverages_pos1\tcontrol_coverages_pos2\tcontrol_contact_distances\tremapping_coverages'
	try:
		if psHash == False: 
			if dups:
				for key in Keys:
					for i in data[key]:
						Str = '%s\t%i\t%s\t%i\t%s\t%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%.2f\t%.5f\n' % (ChrIdxs[key[0]],key[1]*resolution,ChrIdxs[key[2]],key[3]*resolution,i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[9],i[10],i[11])
						print >> f, Str,
			else:
				for key in Keys:
					i = data[key]
					Str = '%s\t%i\t%s\t%i\t%s\t%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%.2f\t%.5f\n' % (ChrIdxs[key[0]],key[1]*resolution,ChrIdxs[key[2]],key[3]*resolution,i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[9],i[10],i[11])
					print >> f, Str,
		else:
			if dups:
				for key in Keys:
					for i in data[key]:
						k = FindKey(psHash,key[3]-key[1])
						Str = '%s\t%i\t%s\t%i\t%s\t%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%.2f\t%.5f\n' % (ChrIdxs[key[0]],key[1]*resolution,ChrIdxs[key[2]],key[3]*resolution,i[0],i[1],psHash[k][i[2]-1],psHash[k][i[3]-1],i[4],i[5],i[6],i[7],i[8],i[9],i[10],i[11])
						print >> f, Str,
			else:
				for key in Keys:
					i = data[key]
					k = FindKey(psHash,key[3]-key[1])
					Str = '%s\t%i\t%s\t%i\t%s\t%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%.2f\t%.5f\n' % (ChrIdxs[key[0]],key[1]*resolution,ChrIdxs[key[2]],key[3]*resolution,i[0],i[1],psHash[k][i[2]-1],psHash[k][i[3]-1],i[4],i[5],i[6],i[7],i[8],i[9],i[10],i[11])
					print >> f, Str,
	except IndexError: print i		
	f.close()
