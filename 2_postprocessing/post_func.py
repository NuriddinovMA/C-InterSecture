import math
import timeit
import sys
import numpy as np
import scipy.stats as st

def boolean(x):
	if x == 'False' or x == 'false': return False
	elif x == 'True' or x == 'true': return True
	else: return x

def HashTry(Hash, key):
	t = 1
	try: Hash[key]
	except KeyError: t = 0
	return t

def colorList(*args):
	try: N = args[0]
	except IndexError: N = 64
	if N < 2: N = 2
	step = int(255/math.ceil(pow(N,1./3)-1))
	color = []
	for i in range(0,256,step):
		for j in range(0,256,step):
			for k in range(0,256,step): color.append( (i/255.,j/255.,k/255.) )
	return color

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

def readBED(file,resolution,type):
	BED = {}
	f1 = open(file, 'r')
	lines = f1.readlines()
	f1.close()
	for i in range(len(lines)-1,0,-1):
		parse = lines[i].split()
		del lines[i]
		try:
			chrName,bin1,bin2 = parse[0], int(parse[1])/resolution, (int(parse[2])-1)/resolution
			try: values = float(parse[3])
			except ValueError: values = parse[3]
			except IndexError: values = '1'
			if bin1 != bin2:
				for i in range(bin1,bin2+1):
					key = chrName,i
					if HashTry(BED, key) == 0: BED[key] = [values,]
					else: BED[key].append(values)
			else:
				key = chrName,bin1
				if HashTry(BED, key) == 0: BED[key] = [values,]
				else: BED[key].append(values)
		except IndexError: break
	Keys = sorted(BED.keys())
	if type == 'values':
		for key in Keys: 
			try: BED[key] = np.nanmean(BED[key])
			except TypeError: BED[key] = len(BED[key])
	else:
		for key in Keys: BED[key] = sorted(set(BED[key]))
	return BED
	
def readContacts(file, Order, resolution, **kwargs):
	try: short = kwargs['short']
	except KeyError: short = False
	try: loci = kwargs['loci']
	except KeyError: loci = False
	start_time = timeit.default_timer()
	f = open(file, 'r')
	lines = f.readlines()
	f.close()
	ln = len(lines)
	if short == False and loci == False:
		Contacts = {}
		for i in range(ln-1,0,-1):
			parse = lines[i].split()
			c0 = int(parse[1])/resolution
			c2 = int(parse[3])/resolution
			if (Order[parse[0]] < Order[parse[2]]): key = parse[0], c0, parse[2], c2
			elif (Order[parse[0]] == Order[parse[2]]) and (c0 <= c2): key = parse[0], c0, parse[2], c2
			else: key = parse[2], c2, parse[0], c0
			c = float(parse[-7])
			q = float(parse[-6])
			dr = float(parse[-5])
			dq =  float(parse[-4])
			l = float(parse[-1])
			if HashTry(Contacts, key) == 0: Contacts[key] = [c,q,dr,dq,l]
			else:
				if l < Contacts[key][-1]: Contacts[key] = [c,q,dr,dq,l]
			if (ln-i) % 1000000 == 0:
				elp = timeit.default_timer() - start_time
				print '\t\tcontact reading progress: %i, time elapsed: %.2f' % (ln-i, elp)
			del lines[i]
	elif short == False and loci == True:
		for i in range(ln-1,0,-1):
			parse = lines[i].split()
			c0 = int(parse[1])/resolution
			c2 = int(parse[3])/resolution
			if (Order[parse[0]] < Order[parse[2]]): key = parse[0], c0, parse[2], c2
			elif (Order[parse[0]] == Order[parse[2]]) and (c0 <= c2): key = parse[0], c0, parse[2], c2
			else: key = parse[2], c2, parse[0], c0
			c = float(parse[-7])
			q = float(parse[-6])
			dr = float(parse[-5])
			dq =  float(parse[-4])
			l = float(parse[-1])
			if HashTry(Contacts, key) == 0: Contacts[key] = [c,q,dr,dq,l]
			else:
				if l < Contacts[key][-1]: Contacts[key] = [c,q,dr,dq,l]
			if (ln-i) % 1000000 == 0:
				elp = timeit.default_timer() - start_time
				print '\t\tcontact reading progress: %i, time elapsed: %.2f' % (ln-i, elp)
			del lines[i]
	else:
		Contacts = []
		for i in range(ln-1,0,-1):
			parse = lines[i].split()
			c0 = int(parse[1])/resolution
			c2 = int(parse[3])/resolution
			if (Order[parse[0]] == Order[parse[2]]): lr = abs(c2-c0)
			else: lr = -1000
			c = float(parse[-7])
			q = float(parse[-6])
			lq = float(parse[-1])
			Contacts.append((c,q,lr,lq))
			if (ln-i) % 1000000 == 0:
				elp = timeit.default_timer() - start_time
				print '\t\tcontact reading progress: %i, time elapsed: %.2f' % (ln-i, elp)
			del lines[i]
	elp = timeit.default_timer() - start_time
	print "\t\ttotal time elapsed:", elp
	return Contacts

def filterLoci(allCon, resolution, loci):
	Keys = allCon.keys()
	ln = len(loci)
	lociCon = [{} for i in range(ln)]
	for key in Keys:
		for l in range(ln):
			if key[0] == loci[l][0] and key[2] == loci[l][0] and loci[l][1]/resolution <= key[1] < key[3] < loci[l][2]/resolution: 
				lociCon[l][key] = allCon[key]
				break
	return lociCon

def randomizeContacts(Contacts,**kwargs):
	shuffle = []
	try: idx = kwargs['idx']
	except KeyError: idx = 1
	try:
		Keys = Contacts.keys()
	except AttributeError:
		randCon = []
		for i in Contacts: shuffle.append(i[idx])
		np.random.shuffle(shuffle)
		for i in range(len(Contacts)):
			D = [Contacts[i][0],Contacts[i][1],Contacts[i][2],Contacts[i][3]]
			D[idx] = shuffle[i]
			randCon.append( D )
	else:
		randCon = {}
		for key in Keys: shuffle.append(Contacts[key][idx])
		np.random.shuffle(shuffle)
		for i in range(len(Keys)):
			D = [Contacts[Keys[i]][0],Contacts[Keys[i]][1],Contacts[Keys[i]][2],Contacts[Keys[i]][3],Contacts[Keys[i]][-1]]
			D[idx] = shuffle[i]
			randCon[Keys[i]] = D
	return randCon

def drawMap(allCon,locus,resolution,Scale1,Scale2,conf):
	Keys = sorted(allCon.keys())
	lb = [ '%s:%ikb' % (locus[0],i*resolution/1000) for i in range(locus[1]/resolution,locus[2]/resolution)]
	Map = [],[],[],[],[],[],lb
	for key in Keys:
		Map[0].append( key[1] - locus[1]/resolution )
		Map[1].append( key[3] - locus[1]/resolution )
		size = (.1*(100 - allCon[key][2]))**2
		size *= Scale1
		color = 100 - allCon[key][0]
		if size > 100: size = 100
		if size < 4: size = 4
		Map[2].append( int(size) )
		Map[3].append( color )
		size = .1*(100-allCon[key][2] - allCon[key][3])
		if size < 2: size = 2
		size = Scale2*size**2
		color = allCon[key][0] - allCon[key][1]
		if size > 100: size = 100
		if color > -20 and color <=0: color = -10
		elif color < 20 and color >0: color = 10
		else: pass
		Map[4].append( int(size) )
		Map[5].append( color )
	return Map

def drawStat(Scale1,Scale2,conf):
	lb = list(range(101))
	Map = [],[],lb
	for i in lb:
		size = (.1*(100 - i))**2
		size *= Scale1
		if size > 100: size = 100
		if size < 4: size = 4
		Map[0].append( int(size) )
		size = .1*(100-i)
		if size < 2: size = 2
		size = Scale2*size**2
		if size > 100: size = 100
		Map[1].append( int(size) )
	return Map

def JuiceboxPre(Contacts,Order,resolution,out):
	Keys = Contacts.keys()
	Keys.sort(key=lambda k:(Order[k[0]],Order[k[2]]))
	f1 = open(out+'.Reference.pre','w')
	f2 = open(out+'.Query.pre','w')
	for key in Keys: 
		print >> f1, '0\t%s\t%i\t0\t1\t%s\t%i\t1\t%f' % (key[0],key[1]*resolution+1,key[2],key[3]*resolution+1, Contacts[key][0])
		print >> f2, '0\t%s\t%i\t0\t1\t%s\t%i\t1\t%f' % (key[0],key[1]*resolution+1,key[2],key[3]*resolution+1, Contacts[key][1])
		del Contacts[key]
	f1.close()
	f2.close()

def trackFeatures(BED):
	Track = {}
	for key in BED:
		feature = BED[key]
		if HashTry(Track, key[0]) == 0: Track[key[0]] = [key[1] for i in range(len(BED[key]))]
		else: Track[key[0]] += [key[1] for i in range(len(BED[key]))]
	for key in Track: Track[key].sort()
	return Track

def trackValues(BED):
	TrackV = {}
	H = {}
	for key in BED:
		values = BED[key]
		if HashTry(H, key[0]) == 0: 
			H[key[0]] = [(key[1],values),]
			TrackV[key[0]] = []
		else: H[key[0]].append((key[1],values))
	for key in H:
		H[key].sort()
		TrackV[key] = [0 for i in range(H[key][-1][0]+1)]
		for i in H[key]: TrackV[key][i[0]] = i[1]
	del H
	return TrackV

def distributionCalc(TrackF,TrackV, frame, repeate):
	Distr = [np.zeros(2*frame+1) for i in range(3)]
	Distr[0] = _distributionCalc(TrackF,TrackV,frame)
	Random = []
	for i in range(repeate):
		TrackFR = randomize(TrackF)
		Random.append(_distributionCalc(TrackFR,TrackV,frame))
	try: Distr[1] = np.nanmean(Random,axis=0)
	except ValueError: Distr[1] = np.zeros(2*frame+1)
	try: Distr[2] = np.nanstd(Random,axis=0,ddof=1)
	except ValueError: Distr[2] = np.zeros(2*frame+1)
	return Distr

def _distributionCalc(TrackF,TrackV,frame):
	Distr = []
	for key in TrackF:
		try: H = TrackV[key]
		except KeyError: pass
		else:
			end = len(H)-1
			mx = max(end,TrackF[key][-1])
			for i in TrackF[key]:
				T = np.zeros(2*frame+1)
				try:
					#if (i-frame) >= 0 and (i+frame+1) <= end: T = H[i-frame:i+frame+1]
					#else: pass
					T = np.zeros(2*frame+1)
					for t in range(i-frame,i+frame+1):
						if t < 0 or t > end: T[t-i+frame] = np.nan
						else: T[t-i+frame] = H[t]
					for t in range(2*frame+1):
						if T[t] == 0: T[t] = np.nan
				except KeyError: pass# print key
				except ValueError: pass# print key,i
				else:
					Distr.append(T)
					#if np.prod(T) != 0: 
					#else: pass
	try: Distr = np.nanmean(Distr,axis=0)
	except ValueError: Distr = np.zeros(2*frame+1)
	return Distr

def randomize(TrackF):
	TrackFR = {}
	for key in TrackF:
		min = TrackF[key][0]
		max = TrackF[key][-1]+1
		size = len(TrackF[key])
		TrackFR[key] = np.random.randint(min,high=max,size=size)
		TrackFR[key].sort()
	return TrackFR

def filterBED(BED1,BED2,threshold,frame):
	BED = {}
	Keys = sorted(BED1)
	for key in Keys:
		if BED1[key] > threshold:
			LocusKeys = [(key[0],i) for i in range(key[1]-frame,key[1]+frame+1)]
			for lk in LocusKeys:
				lk2 = ('chr' + lk[0], lk[1])
				if HashTry(BED, lk) == 0 and HashTry(BED2, lk) == 1: BED[lk] = BED2[lk]
				elif HashTry(BED, lk) == 0 and HashTry(BED2, lk2) == 1: BED[lk] = BED2[lk2]
				else: pass
		else: pass
	return BED

def metricCalc(contacts,resolution,**kwargs):
	Keys = set([])
	metric = []
	try: funckey=kwargs['metric']
	except KeyError: funckey = 'pbad'
	if funckey == 'pearsone': func = _metricCalcPearsone
	elif funckey == 'spearman': func = _metricCalcSpearman
	else: func = _metricCalc
	try: frame=kwargs['frame']
	except KeyError: frame = 8
	try: synblocks=kwargs['synblocks']
	except KeyError: synblocks = False
	try: max_dist=kwargs['max_dist']*1e6/resolution
	except KeyError: max_dist = 100
	try: loci=kwargs['loci']
	except KeyError:
		for key in contacts.keys():
			Keys.add( key[:2] )
			Keys.add( key[2:] )
		Keys = sorted(Keys)
		for key in Keys:
			if synblocks == False or synblocks.has_key(key) == True:
				locus = key[0],key[1],key[1]+1
				locusKeys = [(key[0],i) for i in range(key[1]-frame,key[1]+frame+1)]
				metric.append( func(contacts,locus,locusKeys,max_dist,frame ) )
			else: pass
	else:
		for locus in loci:
			for key in [(locus[0],i) for i in range(locus[1]/resolution,locus[2]/resolution+1)]: locusKeys.add( key )
			locusKeys = sorted(Keys)
			metric.append( func(contacts,locus,locusKeys,max_dist,0) )
	for i in range(len(metric)-1,-1,-1): 
		if metric[i] == None: del metric[i]
	return metric

def _metricCalc(contacts,locus,locusKeys,max_dist,frame):
	I = 0.0
	n = 0
	threshold = 1.5*frame**2
	for j in locusKeys:
		for k in locusKeys:
			if abs(j[1] - k[1]) > max_dist: pass
			else: 
				try: 
					p1 = contacts[j+k][0]
					p2 = contacts[j+k][1]
				except KeyError: pass
				else:
					dp = abs(p1-p2)/100.0
					ds1 = 1.0-abs(p1-50)/50.0
					ds2 = 1.0-abs(p2-50)/50.0
					if ds1 == 0: ds1 = 0.01
					if ds2 == 0: ds2 = 0.01
					if ds1 == 1: ds1 = 0.99
					if ds2 == 1: ds2 = 0.99
					I += -1*dp*np.log10( ds1*ds2 )
					n += 1
	if n > threshold: 
		metric = locus[0],locus[1],locus[2],I/n
		return metric
	else: pass #metric = locus[0],locus[1],locus[2],0

def _metricCalcPearsone(contacts,locus,locusKeys,max_dist,frame):
	I = 0.0
	n = 0
	threshold = 1.5*frame**2
	x = []
	y = []
	for j in locusKeys:
		for k in locusKeys:
			if abs(j[1] - k[1]) > max_dist: pass
			else: 
				try: 
					p1 = contacts[j+k][0]
					p2 = contacts[j+k][1]
				except KeyError: pass
				else:
					x.append(p1)
					y.append(p2)
	if len(x) > threshold: 
		metric = locus[0],locus[1],locus[2],np.corrcoef(x,y)[0,1]
		return metric
	else: pass #metric = locus[0],locus[1],locus[2],0

def _metricCalcSpearman(contacts,locus,locusKeys,max_dist,frame):
	I = 0.0
	n = 0
	threshold = 1.5*frame**2
	x = []
	y = []
	for j in locusKeys:
		for k in locusKeys:
			if abs(j[1] - k[1]) > max_dist: pass
			else: 
				try: 
					p1 = contacts[j+k][0]
					p2 = contacts[j+k][1]
				except KeyError: pass
				else:
					x.append(p1)
					y.append(p2)
	if len(x) > threshold: 
		metric = locus[0],locus[1],locus[2],st.spearmanr(x,y)[0]
		return metric
	else: pass #metric = locus[0],locus[1],locus[2],0
#

#
#multy-species synteny
#
import copy

def readSynBlocks(path,resolution):
	synblocks = []
	f = open(path,'r')
	lines = f.readlines()
	f.close()
	for i in range(len(lines)-1,0,-1):
		parse1 = lines[i].split()
		del lines[i]
		key1 = parse1[0]
		pos1 = int(parse1[1])/resolution, int(parse1[2])/resolution
		if abs(pos1[0] - pos1[1]) > 2: synblocks.append( (key1,pos1[0],pos1[1]) )
	return synblocks

def readSynBlocks2(path,resolution,min_dist):
	synblocks = {}
	f = open(path,'r')
	lines = f.readlines()
	f.close()
	for i in range(len(lines)-1,0,-1):
		parse = lines[i].split()
		del lines[i]
		chrName = parse[0]
		bin1 = int(parse[1])/resolution
		bin2 = int(parse[2])/resolution
		if abs(bin2-bin1) < (2*min_dist + 2): pass
		else:
			for j in range(bin1+min_dist, bin2-min_dist): synblocks[chrName,j] = 1
	return synblocks

def readChrSizes(path,resolution):
	sizes = {}
	f = open(path,'r')
	lines = f.readlines()
	f.close()
	for line in lines:
		parse = line.split()
		sizes[parse[0]] = int(parse[1])/resolution
	del lines
	return sizes

def hashSynBlocks(synblocks,chrSizes):
	hashSynBlocks = {}
	for chrName in chrSizes: hashSynBlocks[chrName] = np.zeros(chrSizes[chrName])
	for sb in synblocks: hashSynBlocks[sb[0]][sb[1]:(sb[2]+1)] = 1
	return hashSynBlocks

def intersectSynBlocks(synblockList, resolution):
	intersect = copy.deepcopy(synblockList[0])
	for sb in synblockList[1:]:
		for key in sb: intersect[key] *= sb[key]
	for key in intersect: intersect[key][1:] -= intersect[key][:-1]
	# print intersect
	for key in intersect.keys():
		c = [],[]
		ln = len(intersect[key])
		for i in range(ln):
			if intersect[key][i] == 1: c[0].append(i)
			if intersect[key][i] == -1: c[1].append(i-1)
		if len(c[1]) < len(c[0]): c[1].append(ln)
		intersect[key] = []
		for i in range(len(c[0])): intersect[key].append( (c[0][i]*resolution,c[1][i]*resolution) )
		if len(intersect[key]) == 0: del intersect[key]
	del c
	return intersect