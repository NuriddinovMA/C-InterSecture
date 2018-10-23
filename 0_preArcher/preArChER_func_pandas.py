import timeit
import sys
try: import numpy as np
except ImportError: print "numpy not found!"
try: import scipy.stats as sc
except ImportError: print "scipy.stats not found!"
try: import pandas
except ImportError: print "pandas not found!"

def HashTry(Hash, key):
	t = 1
	try: Hash[key]
	except KeyError: t = 0
	return t

def unmappedBasesBin(path, resolution, ChrInd, threshold):
	ubh = {}
	f = open(path+'.%i.unmap' % resolution,'r')
	lines = f.readlines()
	f.close()
	for i in range(len(lines)):
		parse = lines[i].split()
		try: 
			unmap = float(parse[3])
			if unmap > 33: ubh[( ChrInd[parse[0]], int(parse[1])/resolution )] = float(parse[3])
			else: pass
		except IndexError: pass
	del lines
	return ubh

def iChromIndexing(path):
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

def iBin2Label(path, ChrInd, resolution):
	binIdxs = [('',0)]
	f = open(path,'r')
	while True:
		parse = f.readline().split()
		try: binIdxs.append( (ChrInd[parse[0]],int(parse[1])/resolution) ) 
		except IndexError: break
	f.close()
	return binIdxs

def iSparseMatrixReader(path, binIdxs, unBinHash, **kwargs):
	try: raw = kwargs['raw']
	except KeyError: raw = True
	try: coverage = kwargs['coverage']
	except KeyError: coverage = 0
	if raw == True: 
		contactHash = _iRawSparceMatrixReader(path,binIdxs,unBinHash,int)
		contactHash = _iContactFiltring(contactHash, unBinHash,coverage)
	else: 
		contactHash = _iNormedSparceMatrixReader(path,binIdxs,float)
		contactHash = _iMatrixNorming(contactHash, unBinHash)
	return contactHash

def _iRawSparceMatrixReader(path,binIdxs,unBinHash,format):
	start_time = timeit.default_timer()
	contactHash = [{},{}]
	f = open(path,'r')
	lines = f.readlines()
	f.close()
	ln = len(lines)
	elp = timeit.default_timer() - start_time
	print '\t\traw sparse matrix lenth: %i contacts, time elapsed: %.2f, memory sized: %.2fMb' % (ln, elp, 1.0*sys.getsizeof(lines)/1024/1024 )
	for i in range(ln-1,-1,-1):
		parse = lines[i].split()
		del lines[i]
		s0 = binIdxs[int(parse[0])]
		s1 = binIdxs[int(parse[1])]
		c = format(parse[2])
		if ( ( (s0[0] == s1[0]) and (abs(s1[1] - s0[1]) > 1) ) or (s0[0] != s1[0]) ):
			if ( HashTry(unBinHash, s0) == 0 ) and ( HashTry(unBinHash, s1) == 0 ): 
				if (s0[0] == s1[0] and s0[1] <= s1[1]) or (s0[0] < s1[0]): 
					if HashTry(contactHash[0], s0+s1) == 0: contactHash[0][s0+s1] = c
					else: contactHash[0][s0+s1] += c
				else:
					if HashTry(contactHash[0], s1+s0) == 0: contactHash[0][s1+s0] = c
					else: contactHash[0][s1+s0] += c
				if HashTry(contactHash[1], s0) == 0: contactHash[1][s0] = c
				else: contactHash[1][s0] += c
				if HashTry(contactHash[1], s1) == 0: contactHash[1][s1] = c
				else: contactHash[1][s1] += c
			else: pass
		else: pass
		if (ln-i) % 10000000 == 0:
			elp = timeit.default_timer() - start_time
			print '\t\traw sparse matrix parsing progress: %i, time elapsed: %.2f, memory sized: %.2fMb' % (ln-i, elp, 1.0*sys.getsizeof(contactHash)/1024/1024 )
	return contactHash

def _iNormedSparceMatrixReader(path,binIdxs,format):
	start_time = timeit.default_timer()
	contactHash = [{},{}]
	f = open(path,'r')
	lines = f.readlines()
	f.close()
	ln = len(lines)
	elp = timeit.default_timer() - start_time
	print '\t\tnormed sparse matrix lenth: %i contacts, time elapsed: %.2f, memory sized: %.2fMb' % (ln, elp, 1.0*sys.getsizeof(lines)/1024/1024 )
	for i in range(ln-1,-1,-1):
		parse = lines[i].split()
		del lines[i]
		s0 = binIdxs[int(parse[0])]
		s1 = binIdxs[int(parse[1])]
		c = format(parse[2])
		if (s0[0] == s1[0] and s0[1] <= s1[1]) or (s0[0] < s1[0]): 
			if HashTry(contactHash[0], s0+s1) == 0: contactHash[0][s0+s1] = c
			else: contactHash[0][s0+s1] += c
		else:
			if HashTry(contactHash[0], s1+s0) == 0: contactHash[0][s1+s0] = c
			else: contactHash[0][s1+s0] += c
		if HashTry(contactHash[1], s0) == 0: contactHash[1][s0] = c
		else: contactHash[1][s0] += c
		if HashTry(contactHash[1], s1) == 0: contactHash[1][s1] = c
		else: contactHash[1][s1] += c
		if (ln-i) % 10000000 == 0:
			elp = timeit.default_timer() - start_time
			print '\t\tnormed sparse matrix parsing progress: %i, time elapsed: %.2f, memory sized: %.2fMb' % (ln-i, elp, 1.0*sys.getsizeof(contactHash)/1024/1024 )
	return contactHash

def _iContactFiltring(contactHash, unBinHash, coverage):
	covList = []
	print '\t\t contact distribution calculation'
	start_time = timeit.default_timer()
	Keys = contactHash[1].keys()
	for key in Keys: covList.append(contactHash[1][key])
	covList.sort()
	p = np.percentile(covList, range(1,101), interpolation='nearest')
	print '\t\t threshold = %i' % coverage
	print '\t\t percentile', p
	del covList
	for key in Keys:
		if contactHash[1][key] <= p[0]: unBinHash[key] = 1
	Keys = contactHash[0].keys()
	for key in Keys:
		if ( HashTry(unBinHash, key[:2]) == 1 ) or ( HashTry(unBinHash, key[2:]) == 1 ) or (contactHash[0][key] == 1):
			c = contactHash[0][key]
			del contactHash[0][key]
			contactHash[1][key[:2]] -= c
			contactHash[1][key[2:]] -= c
			if contactHash[1][key[:2]] < 1: del contactHash[1][key[:2]]
			if contactHash[1][key[2:]] < 1: del contactHash[1][key[2:]]
		else: pass
	del Keys
	elp = timeit.default_timer() - start_time
	print '\t\t contact analyzing full time: %.2f memory sized: %.2fMb' % (elp, 1.0*sys.getsizeof(contactHash[0])/1024/1024 )
	return contactHash

def _iMatrixNorming(contactHash, unBinHash):
	start_time = timeit.default_timer()
	print '\t\t contact distribution calculation'
	l = 0
	Keys = contactHash[0].keys()
	for key in Keys:
		if ( HashTry(unBinHash, key[:2]) == 0 ) and ( HashTry(unBinHash, key[2:]) == 0 ):
			if key[0] == key[2]: l = abs(key[3] - key[1])
			else: l = -1000
			if l == 1 or l == 0: del contactHash[0][key]
			else: contactHash[0][key] = int(contactHash[0][key]/contactHash[1][key[:2]]*1e6),l
		else: del contactHash[0][key]
	del contactHash[1]
	elp = timeit.default_timer() - start_time
	print '\t\t contact analyzing full time: %.2f memory sized: %.2fMb' % (elp, 1.0*sys.getsizeof(contactHash[0])/1024/1024 )
	return contactHash[0]

def iContactBinStat(rawContactHash, normContactHash):
	contactBinList = []
	for key in rawContactHash[0]:
		try:
			contactBinList.append(
				(key[0], key[1], key[2], key[3],
				normContactHash[key][0], rawContactHash[0][key], None, None, rawContactHash[1][key[:2]], rawContactHash[1][key[2:]], normContactHash[key][1])
			)
		except KeyError: pass
	contactBin = pd.DataFrame.from_records(contactBinList, columns=['chr1','pos1','chr2','pos2','contacts','disp','min','max','cov1','cov2','dist'])
	del contactBinList
	contactBin['disp'] = np.sqrt(1./contactBin['disp'])
	return contactBin

def iDistanceHash(pdContactBins,dist_max):
	print '\tDeviation calculate for contactList between %i locus pair' % (len(pdContactBins))
	DistanceHash = {}
	start_time = timeit.default_timer()
	Keys = sorted(pd.unique(pdContactBins['dist']))
	Keys = Keys[1:] + [Keys[0],]
	base_key = Keys[0]
	base_count = pdContactBins[contactBins['dist'] == base_key ]['dist'].count()
	for key in Keys[1:]:
		current_count = pdContactBins[contactBins['dist'] == key]['dist'].count()
		if key >= dist_max or current_count < 0.2*base_count:
			if current_count > base_count/2:
				base_key = key
				DistanceHash[key] = base_key
				sum_count = current_count
			elif sum_count > 0.75*base_count:
				base_key = key
				DistanceHash[key] = base_key
				sum_count = current_count
			else:
				DistanceHash[key] = base_key
				sum_count += current_count
		else: base_key = key
		DistanceHash[key] = base_key
	del Keys
	print '\t', sorted(contactDistance.keys())
	pdContactBins['dist'] = pdContactBins['dist'].applymap(DistanceHash.get)
	elp = timeit.default_timer() - start_time
	print '\tDistance Hashing full time:', elp
	return DistanceHash

def iPercentileStatistics(contactDistanceHash):
	start_time = timeit.default_timer()
	p = range(1,101)
	temp = [i[5] for i in contactDistanceHash]
	try: prcList = np.percentile(temp, p, interpolation='nearest')
	except TypeError:
		#print '\t\tnumpy < 1.9.0',
		prcList = np.percentile(temp, p)
	except NameError:
		#print '\t\tnumpy not found'
		temp.sort()
		ln = len(temp)
		ind = [int(i*ln/100.0) for i in p]
		prcList = [temp[i] for i in ind]
	del temp
	elp = timeit.default_timer() - start_time
	print '\t\ttime to prcList', elp
	return prcList
	
def iPRC(c,prcList):
	prc = [0,0,0]
	try:
		prc[0] = int( sc.percentileofscore(prcList,c[0],kind='rank') )
		prc[1] = int( sc.percentileofscore(prcList,c[0]*(1-c[1]),kind='strict') )
		prc[2] = int( sc.percentileofscore(prcList,c[0]*(1+c[1]),kind='weak') )
	except NameError:
		cmin = c[0]*(1-c[1])
		cmax = c[0]*(1+c[1])
		cmed = [0,0]
		for i in range(len(prcList)):
			if cmin <= prcList[i]:
				prc[1] = i
				break
			else: pass
		for i in range (prc[1],len(prcList)):
			if c[0] <= prcList[i]:
				cmed[0] = i
				break
			else: pass
		else: cmed[0] = len(prcList)
		for i in range (cmed[0],len(prcList)):
			if c[0] > prcList[i]:
				cmed[1] = i
				break
			else: pass
		else: cmed[1] = len(prcList)
		prc[0] = int((cmed[0]+cmed[1])/2)
		for i in range (len(prcList)-1,prc[0],-1):
			if cmax >= prcList[i]:
				prc[2] = i
				break
		else: prc[2] = prc[0]
	return prc

def iPercentilizer(pdContactBins, DistanceHash, prcList):
	prcContactList = []
	for key in DistanceHash['key']: 
		prc = iPRC(pdContactBins[pdContactBins['dist'] == key]['contacts'],prcList)
		pdContactBins['contacts'] = pdContactBins[pdContactBins['dist'] == key]['contacts'].apply(iPRC, args=(prcList))
		pdContactBins['min'] = (pdContactBins[pdContactBins['dist'] == key]['contacts'] * (1. - pdContactBins['disp'])).apply(iPRC, args=(prcList))
		pdContactBins['max'] = (pdContactBins[pdContactBins['dist'] == key]['contacts'] * (1. + pdContactBins['disp'])).apply(iPRC, args=(prcList))
	return prcContactList

def iAbsoluteContacts(pdContactBins):
	absContactList = []
	pdContactBins['min'] = pdContactBins['contacts'] * (1. - pdContactBins['disp'])
	pdContactBins['max'] = pdContactBins['contacts'] * (1. + pdContactBins['disp'])
	return absContactList

def iTotalContactListing(pdContactBins, contactDistanceHash,statistics,resolution,path):
	totalContactList = []
	start_time = timeit.default_timer()
	f = open(path + '.pandas.stat','w')
	Keys = sorted(contactDistanceHash.keys())
	if statistics == 'prc':
		for key in Keys:
			print '\t %i contact transformed for %ikb distance' % (len(contactDistanceHash[key]),key*resolution/1000)
			if key > 0: print >> f, '%ikb\t%i' % (key*resolution/1000, len(contactDistanceHash[key]))
			else: print >> f, 'interchromosome\t%i' % len(contactDistanceHash[key])
			prcList = iPercentileStatistics(contactDistanceHash[key])
			totalContactList += iPercentilizer(contactDistanceHash[key], prcList)
			contactDistanceHash[key] = 0
			del prcList
			elp = timeit.default_timer() - start_time
			print '\t...contact transformed end time: %.2f sec; memory sized: %.2f Mb' % (elp, 1.0*sys.getsizeof(totalContactList)/1024/1024 )
	else:
		for key in Keys:
			print '\t %i contact transformed for %ikb distance' % (len(contactDistanceHash[key]),key*resolution/1000)
			if key > 0: print >> f, '%ikb\t%i' % (key*resolution/1000, len(contactDistanceHash[key]))
			else: print >> f, 'interchromosome\t%i' % len(contactDistanceHash[key])
			totalContactList += iAbsoluteContacts(contactDistanceHash[key])
			contactDistanceHash[key] = 0
			elp = timeit.default_timer() - start_time
			print '\tcontact transformed for %ikb distance for %.2f sec; memory sized: %.2f Mb' % (key*resolution/1000, elp, 1.0*sys.getsizeof(totalContactList)/1024/1024 )
	f.close()
	del Keys
	return totalContactList

def iBinContactWriter(path, contactList, ChrInd):
	print '\tStart contact writing to', path
	start_time = timeit.default_timer()
	contactList.sort()
	elp = timeit.default_timer() - start_time
	print '\tDatabase sorting:',elp
	start_time = timeit.default_timer()
	f = open(path +'.pandas.initialContacts','w')
	print >> f, 'chr1\tpos1\tchr2\tpos2\tcontacts\tmin\tmax\tcov1\tcov2\tdist'
	for i in range(len(contactList)):
		try:
			print >> f, '%s\t%i\t%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i' % (ChrInd[contactList[i][0]],contactList[i][1],ChrInd[contactList[i][2]],contactList[i][3],contactList[i][4],contactList[i][5],contactList[i][6],contactList[i][7],contactList[i][8],contactList[i][9])
		except IndexError:
			print 'index', contactList[i]
		except KeyError:
			print 'key', contactList[i]
	f.close()