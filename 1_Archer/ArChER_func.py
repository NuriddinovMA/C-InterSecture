import math
import timeit
import sys

def HashTry(Hash, key):
	t = 1
	try:
		Hash[key]
	except KeyError:
		t = 0
	return t

def ReadChrIndex(FilePath):
	ChrInd = {}
	File = open(FilePath,'r')
	lines = File.readlines()
	File.close()
	for line in lines:
		parse = line.split('\t')
		ChrInd[parse[0]] = int(parse[1])
		ChrInd[int(parse[1])] = parse[0]
	return ChrInd

def ReadPercentelizedContact(ReadLines): #Reading Contact from database file
	Contact_prc_hash = {}
	start_time = timeit.default_timer()
	for i in range(len(ReadLines)):
		parse = ReadLines[i].split('\t')
		Contact_prc_hash[parse[0], int(parse[1]), parse[2], int(parse[3])] = int(parse[4]),int(parse[5]),int(parse[6]),int(parse[7]),int(parse[8]),int(parse[9])
		if i % 10000000 == 0:
			elp = timeit.default_timer() - start_time
			print '\t\tLocus contact reading:', (i), 'time elapsed:', elp, 'memory sized:', sys.getsizeof(Contact_prc_hash)
	return Contact_prc_hash

def PercentelizedContactStatistic(ReadLines): #Reading Contact from database file
	Stat = [[0 for j in range(10)] for i in range(5)], [[0 for j in range(10)] for i in range(5)], [[0 for j in range(10)] for i in range(5)]
	start_time = timeit.default_timer()
	for i in range(len(ReadLines)):
		parse = ReadLines[i].split('\t')
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
		if i % 10000000 == 0:
			elp = timeit.default_timer() - start_time
			print '\t\tLocus contact reading:', (i), 'time elapsed:', elp, 'memory sized:', sys.getsizeof(Stat)
	for i in range(5):
		for j in range(10):
			Stat[2][i][j] = Stat[1][i][j]/Stat[0][i][j]
	return Stat

def net2mark(name):
	f = open(name, 'r')
	lines = f.readlines()
	f.close()
	M = {}
	L = {}
	for line in lines:
		if line[0] == 'n': 
			key = line.split()[1]
			L[key] = [[],]
		else:
			i = 0
			while line[i] == ' ': i += 1
			ln = (i+1)/2
			if len(L[key]) < ln: L[key].append([line.split()[:7],])
			else: L[key][ln-1].append(line.split()[:7])
	del lines
	chr1 = '-'
	for chrName in L:
		if chr1 != '-':
			end1 = key[2]
			if dir == 1:
				end2 = key[5]
				start2 = end2 - (end1-start1)
			else:
				start2 = key[5]
				end2 = start2 + (end1-start1)
			M[chr1][key].append((start1,end1,start2,end2))
		chr1 = chrName
		M[chr1] = {}
		key = '-'
		for i in L[chr1]:
			for j in i:
				if j[0] == 'fill':
					if key != '-':
						end1 = key[2]
						if dir == 1:
							end2 = key[5]
							start2 = end2 - (end1-start1)
						else:
							start2 = key[5]
							end2 = start2 + (end1-start1)
						M[chr1][key].append((start1,end1,start2,end2))
					start1 = int(j[1])
					end1 = start1+int(j[2])-1
					chr2 = j[3]
					dir = int(j[4]+'1')
					if dir == 1:
						start2 = int(j[5])
						end2 = start2 + int(j[6])-1
					else:
						end2 = int(j[5])
						start2 = end2 - int(j[6])+1
					key = (chr1,start1,end1,chr2,start2,end2,dir)
					M[chr1][key] = []
				else:
					end1 = int(j[1])-1
					if dir == 1:
						end2 = int(j[5])-1
						start2 = end2-(end1-start1)
					else:
						start2 = int(j[5])+1
						end2 = start2+(end1-start1)
					M[chr1][key].append((start1,end1,start2,end2))
					start1 = end1 + int(j[2])+1

	end1 = key[2]
	if dir == 1:
		end2 = key[5]
		start2 = end2 - (end1-start1)
	else:
		start2 = key[5]
		end2 = start2 + (end1-start1)
	M[chr1][key].append((start1,end1,start2,end2))

	del L

	f1 = open(name+'.pre.mark', 'w')
	f2 = open(name+'.2D.ann', 'w')
	print >> f1, 'chr1\tstart1\tend1\tchr2\tstart2\tend2'
	print >> f2, 'chr1\tstart1\tend1\tchr2\tstart2\tend2\tcolor\tcomment'
	Keys = sorted(M.keys())
	for key in Keys:
		syn = sorted(M[key].keys())
		for s in syn:
			if len(s[0]) < 6 and len(s[3]) < 6:
				try: ind = int(s[3][3:])+1
				except ValueError: ind = 0
				if s[2]-s[1] < 40000: 
					c = (s[2]+s[1])/2
					c1,c2 = c - 20000,c + 20000
				else: c1,c2 = s[1],s[2]
				print >> f2, '%s\t%i\t%i\t%s\t%i\t%i\t%s\t%s:%i-%i:%i' % (s[0],c1,c2,s[0],c1,c2,color[ind],s[3],s[4],s[5],s[6])
				for i in M[key][s]: print >> f1, '%s\t%i\t%i\t%s\t%i\t%i' % (s[0],i[0],i[1],s[3],i[2],i[3])
			else: pass
	f1.close()
	f2.close()
	
	M2 = []
	for key in Keys:
		syn = sorted(M[key].keys())
		for s in syn:
			mark = []
			gap1 = 0
			gap2 = 0
			I = len(M[key][s])
			for i in range(I):
				coor = M[key][s][i]
				if len(mark) == 0: mark = [s[0],coor[0],coor[1],s[3],coor[2],coor[3]]
				else:
					gap1 = coor[0]-mark[2]
					gap2 = coor[2]-mark[5]
					if gap1 < 200 and gap2 < 200: mark[2],mark[5] = coor[1],coor[3]
					else: pass
				if gap1 > 200 or gap2 > 200 or i == I-1:
					c = max(mark[5]-mark[4],mark[2]-mark[1])/150
					if c > 1:
						c1,c2 = (mark[2]-mark[1])/c, (mark[5]-mark[4])/c
						for k in range(c-1): M2.append( [s[0],mark[1]+k*c1,mark[1]+k*c1+c1,s[3],mark[4]+k*c2,mark[4]+k*c2+c2] )
						k += 1
						M2.append( [s[0],mark[1]+k*c1,mark[2],s[3],mark[4]+k*c2,mark[5]] )
					else: M2.append( [s[0],mark[1],mark[2],s[3],mark[4],mark[5]] )
					mark = []
				else: pass
	del M
	f1 = open(name+'.mark', 'w')
	print >> f1, 'chr1\tstart1\tend1\tchr2\tstart2\tend2'
	for m in M2:
		if len(m[0]) < 6 and len(m[3]) < 6: print >> f1, '%s\t%i\t%i\t%s\t%i\t%i' % (m[0],m[1],m[2],m[3],m[4],m[5])
		else: pass
	f1.close()
	del M2

def ReadingMarkPoints(lines, Resolution):  #Creating Mark Point List to convert MarkPoint from species to species 
	ObjCoorMP = {},{},{}
	key = ''
	frame = 120000/Resolution
	for line in lines:
		try:
			parse = line.split('\t')
			N1 = parse[0]
			c1 = ( int(parse[1]) + int(parse[2]) ) / 100 * 50
			N2 = parse[3]
			b2 = ( int(parse[4]) + int(parse[5]) ) / (2*Resolution)
	
			if HashTry(ObjCoorMP[1], (N2,b2) ) == 0: ObjCoorMP[1][N2,b2] = set([ (N1,c1) ])
			else: ObjCoorMP[1][N2,b2].add( (N1,c1) )
			if HashTry(ObjCoorMP[2], (N1,c1) ) == 0: ObjCoorMP[2][N1,c1] = set([ (N2,b2) ])
			else: ObjCoorMP[2][N1,c1].add( (N2,b2) )
		except IndexError:  print 'except',line
	
	for i in ObjCoorMP[1]:
		for j in ObjCoorMP[1][i]:
			k = 1.0/len(ObjCoorMP[1][i])
			key = j[0], j[1]/Resolution
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
					else:
						pass
				if p == 0: ObjCoorMP[0][key].append({i:[1,k]})
	for key in ObjCoorMP[0]:
		S = 0.0
		for n in ObjCoorMP[0][key]: 
			for k in n: S += n[k][0]
		for n in ObjCoorMP[0][key]:
			for k in n: n[k] = [ round( n[k][0]/S, 2), round( n[k][1], 2) ]
	sorted(ObjCoorMP[0])
	#return ObjCoorMP
	return ObjCoorMP[0]

def DuplicateContact(Contact_disp_L, ObjCoorMP1, ObjCoorMP2):
	c = [0,0,0,0,0]
	c[4] = ObjCoorMP1[0]*ObjCoorMP2[0]*ObjCoorMP1[1]*ObjCoorMP2[1]
	c[0] = c[4]*Contact_disp_L[0]
	c[1] = c[4]*Contact_disp_L[1]
	c[2] = c[4]*Contact_disp_L[2]
	c[3] = c[4]*Contact_disp_L[-1]
	return c
	
def Label(KeyList, Drive, Resolution):
	lb = ''
	if len(KeyList) < 3 or Drive == 0:
		KeyList = list(KeyList)
		KeyList.sort()	
		for k in KeyList:  
			lb += '%s:%i-' % (k[0],k[1]*Resolution)
	else:
		KeyList = list(KeyList)
		KeyList.sort()
		lb = '%s:%i-' % (KeyList[0][0],KeyList[0][1]*Resolution)
		#for k in KeyList: 
		#	lb += '%s:%ikb-' % (k[0],k[1]*Resolution/1000)
	if lb == '': lb = '-------'
	return lb

	
def DifferContactStat(Contact_disp_0, Contact_disp_1, ObjCoorMP, Model):
	DifferContact = {},{}
	Stat = [ [0 for i in range(101)] for j in range(3)]
	for i in Contact_disp_0:
		key1 = i[:2]
		key2 = i[2:]
		Stat[0][int(Contact_disp_0[i][0])] += 1
		if  (HashTry(ObjCoorMP, key1) == 1) and  (HashTry(ObjCoorMP, key2) == 1):
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

def DifferContact(Contact_disp_0, Contact_disp_1, ObjCoorMP, Resolution, Model, Conf, Equal, Drive):
	DifferContact = [ {} for j in range(7)]
	Stat = [ [0 for i in range(101)] for j in range(8)]
	d = 0
	for i in Contact_disp_0:
		key1 = i[:2]
		key2 = i[2:]
		Stat[0][Contact_disp_0[i][0]/10000] += 1
		if (HashTry(ObjCoorMP, key1) == 1) and  (HashTry(ObjCoorMP, key2) == 1):
			end1 = len(ObjCoorMP[key1])
			end2 = len(ObjCoorMP[key2])
			k = 0
			c = [0,0,0,0,0,set(),set()]
			if end1 == 0: print "!!!", ObjCoorMP[key1]
			if end2 == 0: print "!!!", ObjCoorMP[key2]
			for j1 in range(end1):
				for j2 in range(end2):
					c = [0,0,0,0,0,set(),set()]
					for k1 in ObjCoorMP[key1][j1]:
						for k2 in ObjCoorMP[key2][j2]:
							if HashTry(Contact_disp_1,k1+k2) == 1:
								#print k1+k2, Contact_disp_1[k1+k2], ObjCoorMP[key1][j1][k1], ObjCoorMP[key2][j2][k2]
								dc = DuplicateContact(Contact_disp_1[k1+k2], ObjCoorMP[key1][j1][k1], ObjCoorMP[key2][j2][k2])
								c[0] += dc[0]
								c[1] += dc[1]
								c[2] += dc[2]
								c[3] += dc[3]
								c[4] += dc[4]
								c[5].add(k1)
								c[6].add(k2)
								k += 1
							elif HashTry(Contact_disp_1,k2+k1) == 1:
								#print k2+k1, Contact_disp_1[k2+k1], ObjCoorMP[key1][j1][k1], ObjCoorMP[key2][j2][k2]
								dc = DuplicateContact(Contact_disp_1[k2+k1], ObjCoorMP[key1][j1][k1], ObjCoorMP[key2][j2][k2])
								c[0] += dc[0]
								c[1] += dc[1]
								c[2] += dc[2]
								c[3] += dc[3]
								c[4] += dc[4]
								c[5].add(k2)
								c[6].add(k1)
								k += 1
							else: pass
					if k == 0: 
						Stat[2][Contact_disp_0[i][0]/10000] += 1
						DifferContact[1][i] = [ ( Label(ObjCoorMP[key1][0],0,Resolution)[:-1], Label(ObjCoorMP[key2][0],0,Resolution)[:-1], Contact_disp_0[i][0], Contact_disp_0[i][1], Contact_disp_0[i][2], Contact_disp_0[i][3], Contact_disp_0[i][4], -1 ) ]
					else:
						p = c[4]
						if Model == 1: Norm = 1
						else: Norm = p
						if p != 0: 
							cc = (
								int(round((c[0])/Norm,0)), 
								int(round((c[1])/Norm,0)), 
								int(round((c[2])/Norm,0)),
								int(round((c[3])/Norm,0)),
								set(c[5]),
								set(c[6])
							)
							#if (Contact_disp_0[i][2] < cc[1]) or (Contact_disp_0[i][1] > cc[2]):
							#disp = Contact_disp_0[i][2] - Contact_disp_0[i][1] + cc[2]  - cc[1]
							disp1 = max(( Contact_disp_0[i][2]-Contact_disp_0[i][0]),(Contact_disp_0[i][0]-Contact_disp_0[i][1])) 
							disp2 = max((cc[2]-cc[0]),(cc[0]-cc[1]))
							disp = disp1+disp2+1
							if (disp1+disp2) < 0: print '\t\tError', i, k1, k2
							dfr = math.fabs(Contact_disp_0[i][0] - cc[0])
							conf = round(1.0*dfr/disp, 7)
							Stat[4][Contact_disp_0[i][0]/10000] += 1

							if HashTry(DifferContact[3],i) == 0: DifferContact[3][i] = [ (Label(cc[4],0,Resolution)[:-1], Label(cc[5],0,Resolution)[:-1], Contact_disp_0[i][0], cc[0], disp1, disp2, Contact_disp_0[i][3], Contact_disp_0[i][4], cc[3]) ]
							else: DifferContact[3][i].append( (Label(cc[4],0,Resolution)[:-1], Label(cc[5],0,Resolution)[:-1], Contact_disp_0[i][0], cc[0], disp1, disp2, Contact_disp_0[i][3], Contact_disp_0[i][4], cc[3]) )

							if Drive == 0:
								if conf > Conf:
									if HashTry(DifferContact[5],i) == 0: DifferContact[5][i] = [ (Label(cc[4],0,Resolution)[:-1], Label(cc[5],0,Resolution)[:-1], Contact_disp_0[i][0], cc[0], disp1, disp2, Contact_disp_0[i][3], Contact_disp_0[i][4], cc[3]) ]
									else: DifferContact[5][i].append( (Label(cc[4],0,Resolution)[:-1], Label(cc[5],0,Resolution)[:-1], Contact_disp_0[i][0], cc[0], disp1, disp2, Contact_disp_0[i][3], Contact_disp_0[i][4], cc[3]) )
								elif (disp <= Equal) and (dfr < disp):
									if HashTry(DifferContact[6],i) == 0: DifferContact[6][i] = [ (Label(cc[4],0,Resolution)[:-1], Label(cc[5],0,Resolution)[:-1], Contact_disp_0[i][0], cc[0], disp1, disp2, Contact_disp_0[i][3], Contact_disp_0[i][4], cc[3]) ]
									else: DifferContact[6][i].append( (Label(cc[4],0,Resolution)[:-1], Label(cc[5],0,Resolution)[:-1], Contact_disp_0[i][0], cc[0], disp1, disp2, Contact_disp_0[i][3], Contact_disp_0[i][4], cc[3]) )
								else: pass
							else:
								ccl = list(cc[4]),list(cc[5])
								ccl[0].sort()
								ccl[1].sort()
								k = 0
								dlen = math.fabs( math.fabs(i[1]-i[3]) - math.fabs(ccl[0][0][1] - ccl[1][-1][1]) )
								if (conf > Conf) and (i[0] == i[2]) and (ccl[0][0][0] == ccl[1][-1][0]) and (dlen < (5000000/Resolution)):
									if HashTry(DifferContact[5],i) == 0: DifferContact[5][i] = [ (Label(cc[4],0,Resolution)[:-1], Label(cc[5],0,Resolution)[:-1], Contact_disp_0[i][0], cc[0], disp1, disp2, Contact_disp_0[i][3], Contact_disp_0[i][4], cc[3]) ]
									else: DifferContact[5][i].append( (Label(cc[4],0,Resolution)[:-1], Label(cc[5],0,Resolution)[:-1], Contact_disp_0[i][0], cc[0], disp1, disp2, Contact_disp_0[i][3], Contact_disp_0[i][4], cc[3]) )
									Stat[6][Contact_disp_0[i][0]/10000] += 1
									k += 1
								elif (disp <= Equal) and (dfr < disp) and (i[0] == i[2]) and (ccl[0][0][0] == ccl[1][-1][0]) and (dlen < (5000000/Resolution)):
									if HashTry(DifferContact[6],i) == 0: DifferContact[6][i] = [ (Label(cc[4],0,Resolution)[:-1], Label(cc[5],0,Resolution)[:-1], Contact_disp_0[i][0], cc[0], disp1, disp2, Contact_disp_0[i][3], Contact_disp_0[i][4], cc[3]) ]
									else: DifferContact[6][i].append( (Label(cc[4],0,Resolution)[:-1], Label(cc[5],0,Resolution)[:-1], Contact_disp_0[i][0], cc[0], disp1, disp2, Contact_disp_0[i][3], Contact_disp_0[i][4], cc[3]) )
									Stat[7][Contact_disp_0[i][0]/10000] += 1
									k += 1
								else: pass  
								if k == 0:
									Stat[5][Contact_disp_0[i][0]/10000] += 1
									DifferContact[4][i] =  [(Label(cc[4],0,Resolution)[:-1], Label(cc[5],0,Resolution)[:-1], Contact_disp_0[i][0], cc[0], disp1, disp2, Contact_disp_0[i][3], Contact_disp_0[i][4], cc[3])]
								else: pass
						else:
							Stat[3][Contact_disp_0[i][0]/10000] += 1
							if HashTry(DifferContact[2],i) == 0: DifferContact[2][i] = [ ( Label(ObjCoorMP[key1][0],0,Resolution)[:-1], Label(ObjCoorMP[key2][0],0,Resolution)[:-1], Contact_disp_0[i][0], -1, -1, Contact_disp_0[i][3], Contact_disp_0[i][4], -1 ) ]
							else: DifferContact[2][i].append( ( Label(ObjCoorMP[key1][0],0,Resolution)[:-1], Label(ObjCoorMP[key2][0],0,Resolution)[:-1], Contact_disp_0[i][0], -1, -1, Contact_disp_0[i][3], Contact_disp_0[i][4], -1 ) )
		else: 
			Stat[1][Contact_disp_0[i][0]/10000] += 1
			DifferContact[0][i] =  [('unmapped', 'unmapped', Contact_disp_0[i][0], Contact_disp_0[i][1], Contact_disp_0[i][2], Contact_disp_0[i][3], Contact_disp_0[i][4], -1)]
	for i in range(8): Stat[i][0] = sum(Stat[i])
	print d
	print '\t\t Initial pair locus:', Stat[0][0], 'for 10%-groups:', Stat[0]
	print '\t\t Unmapped pair locus:', Stat[1][0], 'for 10%-groups:', Stat[1]
	print '\t\t Uncontacted pair locus:', Stat[2][0], 'for 10%-groups:', Stat[2]
	print '\t\t Map-covering pair locus:', Stat[3][0], 'for 10%-groups:', Stat[3]
	print '\t\t All filtred pair locus:', Stat[4][0], 'for 10%-groups:', Stat[4]
	print '\t\t Nonsingnificant pair locus:', Stat[5][0], 'for 10%-groups:', Stat[5]
	print '\t\t Differ pair locus:', Stat[6][0], 'for 10%-groups:', Stat[6]
	print '\t\t Equal pair locus:', Stat[7][0], 'for 10%-groups:', Stat[7]
	return DifferContact,Stat

def DifferContactShort(Contact_disp_0, Contact_disp_1, ObjCoorMP, Resolution, Model, Conf, Equal, Drive):
	DifferContact = {}

	for i in Contact_disp_0:
		key1 = i[:2]
		key2 = i[2:]
		if (HashTry(ObjCoorMP, key1) == 1) and  (HashTry(ObjCoorMP, key2) == 1):
			end1 = len(ObjCoorMP[key1])
			end2 = len(ObjCoorMP[key2])
			k = 0
			c = [0,0,0,0,0,set(),set()]
			if end1 == 0: print "!!!", ObjCoorMP[key1]
			if end2 == 0: print "!!!", ObjCoorMP[key2]
			for j1 in range(end1):
				for j2 in range(end2):
					c = [0,0,0,0,0,set(),set()]
					for k1 in ObjCoorMP[key1][j1]:
						for k2 in ObjCoorMP[key2][j2]:
							if HashTry(Contact_disp_1,k1+k2) == 1:
								#print k1+k2, Contact_disp_1[k1+k2], ObjCoorMP[key1][j1][k1], ObjCoorMP[key2][j2][k2]
								dc = DuplicateContact(Contact_disp_1[k1+k2], ObjCoorMP[key1][j1][k1], ObjCoorMP[key2][j2][k2])
								c[0] += dc[0]
								c[1] += dc[1]
								c[2] += dc[2]
								c[3] += dc[3]
								c[4] += dc[4]
								c[5].add(k1)
								c[6].add(k2)
								k += 1
							elif HashTry(Contact_disp_1,k2+k1) == 1:
								#print k2+k1, Contact_disp_1[k2+k1], ObjCoorMP[key1][j1][k1], ObjCoorMP[key2][j2][k2]
								dc = DuplicateContact(Contact_disp_1[k2+k1], ObjCoorMP[key1][j1][k1], ObjCoorMP[key2][j2][k2])
								c[0] += dc[0]
								c[1] += dc[1]
								c[2] += dc[2]
								c[3] += dc[3]
								c[4] += dc[4]
								c[5].add(k2)
								c[6].add(k1)
								k += 1
							else: pass
					if k == 0: pass
					else:
						p = c[4]
						if Model == 1: Norm = 1
						else: Norm = p
						if p != 0: 
							cc = (
								int(round((c[0])/Norm,0)), 
								int(round((c[1])/Norm,0)), 
								int(round((c[2])/Norm,0)),
								int(round((c[3])/Norm,0)),
								set(c[5]),
								set(c[6])
							)
							disp1 = max(( Contact_disp_0[i][2]-Contact_disp_0[i][0]),(Contact_disp_0[i][0]-Contact_disp_0[i][1])) 
							disp2 = max((cc[2]-cc[0]),(cc[0]-cc[1]))
							disp = disp1+disp2+1
							if (disp1+disp2) < 0: print '\t\tError', i, k1, k2
							dfr = math.fabs(Contact_disp_0[i][0] - cc[0])
							conf = round(1.0*dfr/disp, 7)


							if HashTry(DifferContact,i) == 0: DifferContact[i] = [ (Label(cc[4],0,Resolution)[:-1], Label(cc[5],0,Resolution)[:-1], Contact_disp_0[i][0], cc[0], disp1, disp2, Contact_disp_0[i][3], Contact_disp_0[i][4], cc[3]) ]
							else: DifferContact[i].append( (Label(cc[4],0,Resolution)[:-1], Label(cc[5],0,Resolution)[:-1], Contact_disp_0[i][0], cc[0], disp1, disp2, Contact_disp_0[i][3], Contact_disp_0[i][4], cc[3]) )

						else: pass
							
		else: pass

	return DifferContact

def PrintDifferContact(Inf, Resolution, out):
	f = open(out, 'w')
	KeysList = Inf.keys()
	KeysList.sort()
	print >> f, 'chr1\tpos1\tchr2\tpos2\tremap1\tremap2\tobs\texp\tdisp_obs\tdisp_exp\tcvr1\tcvr2\texp_dist'
	for Keys in KeysList:
		for i in Inf[Keys]:
			try:
				Str = '%s\t%i\t%s\t%i\t%s\t%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n' % (Keys[0],Keys[1]*Resolution,Keys[2],Keys[3]*Resolution,i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8])
			except:	
				Str = '%s\t%i\t%s\t%i\t%s\t%s\t%i\t%i\t%i\t%i\t%i\n' % (Keys[0],Keys[1]*Resolution,Keys[2],Keys[3]*Resolution,i[0],i[1],i[2],i[3],i[4],i[5],i[6])
			print >> f, Str,
	f.close()
	return out