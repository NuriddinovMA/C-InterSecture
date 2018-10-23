import math
import timeit
import sys
import numpy as np

def HashTry(Hash, key):
	t = 1
	try: Hash[key]
	except KeyError: t = 0
	return t

def chrSortOrder(chrSizes):
	order = {}
	file = open(chrSizes,'r')
	lines = file.readlines()
	file.close()
	for l in range(len(lines)): order[lines[l].split()[0]] = l
	return order

def readBED(file,resolution,type):
	BED = {}
	f1 = open(file, 'r')
	line = f1.readline()
	while True:
		line = f1.readline()
		parse = line.split()
		try: 
			chrName,bin1,bin2 = parse[0], int(parse[1])/resolution, (int(parse[2])-1)/resolution
			try: values = float(parse[3])
			except ValueError: values = parse[3]
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
	f1.close()
	Keys = sorted(BED.keys())
	if type == 'values':
		for key in Keys: 
			try: BED[key] = np.nanmean(BED[key])
			except TypeError: BED[key] = len(BED[key])
	else:
		for key in Keys: BED[key] = sorted(set(BED[key]))
	return BED

def readContacts(file, resolution):
	start_time = timeit.default_timer()
	Contacts = {}
	f1 = open(file, 'r')
	lines = f1.readlines()
	f1.close()
	for i in range(len(lines)-1,0,-1):
		parse = lines[i].split()
		key = parse[0], int(parse[1])/resolution, parse[2], int(parse[3])/resolution
		d = float(parse[-5]) + float(parse[-4]) + 1
		c = float(parse[-7])
		f = float(parse[-6])
		v = 1.0*abs(c-f)/d
		l = float(parse[-1])
		if HashTry(Contacts, key) == 0: Contacts[key] = [c,f,d,v,l]
		else:
			if l < Contacts[key][-1]: Contacts[key] = [c,f,d,v,l]
	elp = timeit.default_timer() - start_time
	print "time elapsed 1:", elp
	del lines
	elp = timeit.default_timer() - start_time
	print "time elapsed 2:", elp
	return Contacts

def readContacts2(file, Order, resolution):
	start_time = timeit.default_timer()
	Contacts = {}
	f1 = open(file, 'r')
	line = f1.readline()
	while True:
		line = f1.readline()
		parse = line.split()
		try:
			if Order[parse[0]] < Order[parse[2]]: key = parse[0], int(parse[1])/resolution, parse[2], int(parse[3])/resolution
			elif parse[0] == parse[2]:
				c1 = int(parse[3])/resolution
				c2 = int(parse[1])/resolution
				if c1 < c2: key = parse[0], c1, parse[2], c2
				else: key = parse[0], c2, parse[2], c1
			else: key = parse[2], int(parse[3])/resolution, parse[0], int(parse[1])/resolution
		except IndexError: break
		d = float(parse[-5]) + float(parse[-4]) + 1
		c = float(parse[-7])
		f = float(parse[-6])
		v = 1.0*abs(c-f)/d
		l = float(parse[-1])
		if HashTry(Contacts, key) == 0: Contacts[key] = [c,f,d,v,l]
		else:
			if l < Contacts[key][-1]: Contacts[key] = [c,f,d,v,l]
	f1.close()
	elp = timeit.default_timer() - start_time
	print "\t\ttime elapsed:", elp
	return Contacts

def JuiceboxPre(Contacts,Order,resolution,out):
	Keys = Contacts.keys()
	Keys.sort(key=lambda k:(Order[k[0]],Order[k[2]],k[1],k[3]))
	f = open(out+'.Reference.pre','w')
	for key in Keys: print >> f, '0\t%s\t%i\t0\t1\t%s\t%i\t1\t%f' % (key[0],key[1]*resolution+1,key[2],key[3]*resolution+1, Contacts[key][0])
	f.close()
	f = open(out+'.Query.pre','w')
	for key in Keys: print >> f, '0\t%s\t%i\t0\t1\t%s\t%i\t1\t%f' % (key[0],key[1]*resolution+1,key[2],key[3]*resolution+1, Contacts[key][1])
	f.close()

def trackValues(BED):
	Track = {}
	H = {}
	for key in BED:
		values = BED[key]
		if HashTry(H, key[0]) == 0: 
			H[key[0]] = [(key[1],values),]
			Track[key[0]] = []
		else: H[key[0]].append((key[1],values))
	for key in H: 
		H[key].sort()
		Track[key] = [0 for i in range(H[key][-1][0]+1)]
		for i in H[key]: Track[key][i[0]] = i[1]
	del H
	return Track

def trackFeatures(BED):
	Track = {}
	for key in BED:
		feature = BED[key]
		if HashTry(Track, key[0]) == 0: Track[key[0]] = [key[1],]
		else: Track[key[0]].append(key[1])
	for key in Track: Track[key].sort()
	return Track

def distributionCalc(TrackF,TrackV,frame,repeate):
	Distr = [np.zeros(2*frame+1) for i in range(3)]
	Distr[0] = _distributionCalc(TrackF,TrackV,frame)
	Random = []
	for i in range(repeate):
		TrackFR = randomize(TrackF)
		Random.append(_distributionCalc(TrackFR,TrackV,frame))
	Distr[1] = np.mean(Random,axis=0)
	Distr[2] = np.std(Random,axis=0,ddof=1)
	return Distr

def _distributionCalc(TrackF,TrackV,frame):
	Distr = np.zeros(2*frame+1)
	n = 0
	#mean = 0
	for key in TrackF:
		#mean += np.mean(TrackV[key])
		for i in TrackF[key]:
			try: 
				T = TrackV[key][(i-frame):(i+frame+1)]
				h = 0
				for t in T: 
					if t == 0: h += 1
				if h < 5:
					Distr += T
					n += 1
				else: pass
			except KeyError: pass# print key
			except ValueError: pass# print key,i
	Distr = Distr / n
	#mean = 1.0*mean/len(TrackF)
	return Distr

def randomize(TrackF):
	TrackFR = {}
	for key in TrackF:
		min = TrackF[key][0]
		max = TrackF[key][-1]-min+1
		size = len(TrackF[key])
		TrackFR[key] = min+np.random.choice(max,size=size,replace=False)
		TrackFR[key].sort()
		#if key == 'chr19': print TrackFR[key]
	return TrackFR

def filterBED(BED1,BED2,threshold,frame):
	BED = {}
	Keys = sorted(BED1)
	for key in Keys:
		if BED1[key] > threshold:
			LocusKeys = [(key[0],i) for i in range(key[1]-frame,key[1]+frame+1)]
			for lk in LocusKeys:
				if HashTry(BED, lk) == 0 and HashTry(BED2, lk) == 1: BED[lk] = BED2[lk]
				else: pass
		else: pass
	return BED

def metricCalc(contacts,resolution,**kwargs):
	Keys = set([])
	metric = []
	try: frame=kwargs['frame']
	except KeyError: frame = 8
	try: max_dist=kwargs['max_dist']*1e6/resolution
	except KeyError: max_dist = 100
	try: loci=kwargs['loci']
	except KeyError:
		for key in contacts.keys():
			Keys.add( key[:2] )
			Keys.add( key[2:] )
		Keys = sorted(Keys)
		for key in Keys:
			locus = key[0],key[1],key[1]+1
			locusKeys = [(key[0],i) for i in range(key[1]-frame,key[1]+frame+1)]
			metric.append( _metricCalc(contacts,locus,locusKeys,frame,max_dist) )
	else:
		for locus in loci:
			for key in [(locus[0],i) for i in range(locus[1]/resolution,locus[2]/resolution+1)]: locusKeys.add( key )
			locusKeys = sorted(Keys)
			metric.append( _metricCalc(contacts,locus,locusKeys,frame,max_dist) )
	for i in range(len(metric)-1,-1,-1): 
		if metric[i] == None: del metric[i]
	return metric

def _metricCalc(contacts,locus,locusKeys,frame,max_dist):
	I = 0.0
	n = 0
	threshold = (frame-3)*frame/4
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