import math
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

def iReadPercentelizedContact(path,ChrIdxs): #Reading Contact from database file
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

def netParser(name):
	f = open(name, 'r')
	lines = f.readlines()
	f.close()
	parsedNet = {}
	for line in lines:
		if line[0] == 'n': 
			key = line.split()[1]
			line = 0
			parsedNet[key] = [[],]
		else:
			k = 0
			while line[k] == ' ': k += 1
			ln = (k+1)/2
			parse = line.split()[:7]
			data = parse[0],int(parse[1]),int(parse[2]),parse[3],int(parse[5]),int(parse[6]),int(parse[4]+'1')
			if len(parsedNet[key]) < ln: parsedNet[key].append([data,])
			else: parsedNet[key][ln-1].append(data)
			line = 0
	del lines
	return parsedNet

def net2pre(parsedNet,name):
	preMark = {}
	for chr1 in parsedNet:
		preMark[chr1] = {}
		for i in parsedNet[chr1]:
			for j in i:
				start1 = j[1]
				end1 = start1+j[2]-1
				chr2 = j[3]
				dir = j[6]
				if j[0] == 'fill':
					start2 = j[4]
					end2 = j[4] + j[5] - 1
					if dir == 1: pass
					else: start2,end2 = end2,start2
					key = (chr1,start1,end1,chr2,start2,end2,dir)
					preMark[chr1][key] = [(start1,start1,start2,start2),]
				else:
					start2 = j[4] - 1*dir
					end2 = j[4] + j[5]*dir
					preMark[chr1][key].append((start1,end1,start2,end2))
	for chrName in preMark:
		for key in preMark[chrName]:
			m = preMark[chrName][key]
			temp = []
			if len(m) == 1: preMark[chrName][key] = [(key[1],key[2],key[4],key[5])]
			else:
				for i in range(1,len(m)): temp.append( (m[i-1][1],m[i][0],m[i-1][3],m[i][2]) )
				temp.append( (m[i][1],key[2],m[i][3],key[5]) )
				preMark[chrName][key] = temp
	color = colorList()
	f1 = open(name+'.pre.mark', 'w')
	f2 = open(name+'.2D.ann', 'w')
	print >> f1, 'chr1\tstart1\tend1\tchr2\tstart2\tend2'
	print >> f2, 'chr1\tstart1\tend1\tchr2\tstart2\tend2\tcolor\tcomment'
	chrNames = sorted(preMark.keys())
	for chrName in chrNames:
		syn = sorted(preMark[chrName].keys())
		for key in syn:
			if len(key[0]) < 6 and len(key[3]) < 6:
				try: ind = int(key[3][3:])+1
				except ValueError: ind = 0
				if key[2]-key[1] < 5000: 
					c = (key[2]+key[1])/2
					c1,c2 = c - 2500,c + 2500
				else: c1,c2 = key[1],key[2]
				#print ind,color[ind]
				print >> f2, '%s\t%i\t%i\t%s\t%i\t%i\t%s\t%s:%i-%i:%i' % (key[0],c1,c2,key[0],c1,c2,color[ind],key[3],key[4],key[5],key[6])
				for i in preMark[chrName][key]: print >> f1, '%s\t%i\t%i\t%s\t%i\t%i' % (key[0],i[0],i[1],key[3],i[2],i[3])
			else: pass
	f1.close()
	f2.close()
	return preMark

def pre2mark(preMark,name):
	markPoints = []
	chrNames = sorted(preMark.keys())
	for chrName in chrNames:
		syn = sorted(preMark[chrName].keys())
		for s in syn:
			ln = len(preMark[chrName][s])
			gap = 0,0
			if ln == 1:
				coor = preMark[chrName][s][0]
				mark = [s[0],coor[0],coor[1],s[3],coor[2],coor[3]]
				l1,l2 = mark[2]-mark[1],mark[5]-mark[4]
				c = min( l1, abs(l2) )/200
				if c > 0: 
					c1,c2 = 1.*l1/(c+1),1.*l2/(c+1)
					for k in range(c+1): markPoints.append( [mark[0], int(mark[1] + c1*k), int(mark[1] + c1*(k+1)), mark[3], int(mark[4] + c2*k), int(mark[4] + c2*(k+1)) ] )
				else: markPoints.append(mark)
			else:
				for i in range(ln):
					coor = preMark[chrName][s][i]
					if i == 0: mark = [s[0],coor[0],coor[1],s[3],coor[2],coor[3]]
					elif i == ln-1:
						mark[2],mark[5] = coor[1],coor[3]
						l1,l2 = mark[2]-mark[1],mark[5]-mark[4]
						c = min( l1, abs(l2) )/200
						if c > 0: 
							c1,c2 = 1.*l1/(c+1),1.*l2/(c+1)
							for k in range(c+1): markPoints.append( [mark[0], int(mark[1] + c1*k), int(mark[1] + c1*(k+1)), mark[3], int(mark[4] + c2*k), int(mark[4] + c2*(k+1)) ] )
						else: markPoints.append(mark)
					else:
						gap = coor[0]-mark[2],abs(coor[2]-mark[5])
						if gap[0] < 200 and gap[1] < 200: mark[2],mark[5] = coor[1],coor[3]
						else:
							l1,l2 = mark[2]-mark[1],mark[5]-mark[4]
							c = min( l1, abs(l2) )/200
							if c > 0:
								c1,c2 = 1.*l1/(c+1),1.*l2/(c+1)
								for k in range(c+1): markPoints.append( [mark[0], int(mark[1] + c1*k), int(mark[1] + c1*(k+1)), mark[3], int(mark[4] + c2*k), int(mark[4] + c2*(k+1)) ] )
							else: markPoints.append(mark)
							mark = [ s[0],coor[0],coor[1],s[3],coor[2],coor[3] ]
	f1 = open(name+'.mark', 'w')
	print >> f1, 'chr_sp1\tstart1\tend1\tchr_sp2\tstart2\tend2'
	for m in markPoints:
		if len(m[0]) < 6 and len(m[3]) < 6 and (m[2]-m[1] > 15) and (abs(m[5] - m[4]) > 15): print >> f1, '%s\t%i\t%i\t%s\t%i\t%i' % (m[0],m[1],m[2],m[3],m[4],m[5])
		else: pass
	f1.close()
	del markPoints

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
	else:
		if (writed[-1] < dupled[-1]) or (dupled[-1] < 0): results = True
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
	Statistic = ['all','remappable','processed','uniq','duplicated','dropped'],[0,0,0,set([]),[],0],[set([]),set([]),set([]),set([]),set([]),set([])]
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
							if (k1[0] != k2[0]) and inter == False: pass
							elif HashTry(Contact_disp_1,k1+k2) == 1:
								dc = iDuplicateContact(Contact_disp_1[k1+k2], ObjCoorMP[key1][j1][k1], ObjCoorMP[key2][j2][k2])
								c[0] += dc[0]
								c[1] += dc[1]
								c[2] += dc[2]
								c[3] += dc[3]
								c[4] += dc[4]
								c[-4] += dc[-4]
								c[-3] += dc[-3]
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
								c[-4] += dc[-4]
								c[-3] += dc[-3]
								c[-2].add(k2)
								c[-1].add(k1)
								k += 1
							else: Statistic[2][5] |= set([key1,key2,])
					if k == 0: pass
					else:
						Statistic[1][2] += 1
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
								set(c[-2]),
								set(c[-1])
							)
							disp1 = max(( Contact_disp_0[i][2]-Contact_disp_0[i][0]),(Contact_disp_0[i][0]-Contact_disp_0[i][1])) 
							disp2 = max((cc[2]-cc[0]),(cc[0]-cc[1]))
							to_write = (iLabel(cc[-2],resolution,ChrIdxs2)[:-1], iLabel(cc[-1],resolution,ChrIdxs2)[:-1], Contact_disp_0[i][0], cc[0], disp1, disp2, Contact_disp_0[i][3], Contact_disp_0[i][4], cc[3], cc[4], cc[-3])
							if HashTry(DifferContact,i) == 0: 
								DifferContact[i] = to_write
								Statistic[1][3] |= set([i])
							else:
								Statistic[1][4].append(i)
								#Statistic[2][4] += [key1,key2,] 
								#if DifferContact[i][-1] < cc[-3] or cc[-3] <= 0 : pass
								if iChouseBest(DifferContact[i],to_write,criteria) == True: pass
								else: DifferContact[i] = to_write
						else: pass
		else: pass
	Statistic[1][5] = len(Statistic[1][4])
	Statistic[2][5] = Statistic[2][5] - Statistic[2][2]
	Statistic[2][3] = Statistic[2][3] - Statistic[2][5]
	Statistic[2][4] = Statistic[2][4] - Statistic[2][5]
	N = len( Statistic[1][3] & set(Statistic[1][4]) )
	Statistic[1][3] = len(Statistic[1][3]) - N
	Statistic[1][4] = Statistic[1][5] + N
	f = open(stat_out+'.stat','w')
	for i in range(6): 
		Statistic[2][i] = len(Statistic[2][i])
		print >> f, Statistic[0][i],Statistic[1][i],Statistic[2][i]
	del Statistic
	return DifferContact

def iPrintDifferContact(data, resolution, ChrIdxs, out):
	f = open(out, 'w')
	Keys = sorted(data.keys(),key=lambda x: (x[0],x[2],x[1],x[3]))
	print >> f, 'chr1_observed\tpos1_observed\tchr2_observed\tpos2_observed\tremap1_control\tremap2_control\tobserved_contacts\tcontrol_contacts\tobserved_deviations\tcontrol_deviations\tobserved_coverages_pos1\tobserved_coverages_pos2\tcontrol_coverages_pos1\tcontrol_coverages_pos2\tcontrol_contact_distances'
	for key in Keys:
		i = data[key]
		#for i in data[key]:
		try:
			Str = '%s\t%i\t%s\t%i\t%s\t%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%.2f\n' % (ChrIdxs[key[0]],key[1]*resolution,ChrIdxs[key[2]],key[3]*resolution,i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[9],i[10])
		except:
			Str = '%s\t%i\t%s\t%i\t%s\t%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n' % (ChrIdxs[key[0]],key[1]*resolution,ChrIdxs[key[2]],key[3]*resolution,i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8])
		print >> f, Str,
	f.close()
