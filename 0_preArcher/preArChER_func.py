import math
import timeit
import sys
import numpy as np

def HashTry(Hash, key):
	t = 1
	try: Hash[key]
	except KeyError: t = 0
	return t

def iReadBin2Label(ReadLines, ChrInd, Bin_size):
	iRB2L_list = [('',0)]
	for line in ReadLines:
		parse = line.split()
		iRB2L_list.append( (ChrInd[parse[0]],int(parse[1])/Bin_size) ) 
	return iRB2L_list

def iRawMatrixFilter(ReadLines, iRB, th, unBinHash):
	filtredContact = {}
	Cov = []
	print '\tStart filtring. Contact analysis...'
	start_time = timeit.default_timer()
	for i in range(len(ReadLines)):
		parse = ReadLines[i].split()
		s0 = iRB[int(parse[0])]
		s1 = iRB[int(parse[1])]
		c = int(parse[2])

		if ( ( (s0[0] == s1[0]) and (math.fabs(s1[1] - s0[1]) > 1) ) or (s0[0] != s1[0]) ):
			if ( HashTry(unBinHash, s0) == 0 ) and ( HashTry(unBinHash, s1) == 0 ):
				if HashTry(filtredContact, s0) == 0: filtredContact[s0] = c
				else: filtredContact[s0] += c
				if HashTry(filtredContact, s1) == 0: filtredContact[s1] = c
				else: filtredContact[s1] += c
	
		if i % 10000000 == 0:
			elp = timeit.default_timer() - start_time
			print '\t\tContact filtring progress: %i, time elapsed: %.2f, memory sized: %.2f' % ( i, elp, sys.getsizeof(filtredContact) )
	
	Keys = filtredContact.keys()
	for key in Keys:
		Cov.append(filtredContact[key])
	Cov.sort()
	thold = Cov[int(th*len(Cov)/100.0)]
	print '\t\t threshold = %i' % thold
	print '\t\t percentile', np.percentile(Cov, range(1,101))
	del Cov
	Keys = filtredContact.keys()
	for key in Keys:
		if filtredContact[key] > thold: del filtredContact[key]

	elp = timeit.default_timer() - start_time
	print '\tContact analyzing full time: %.2f memory sized: %.2f' % (elp, sys.getsizeof(filtredContact) )
	return filtredContact

def iRawMatrixRead(ReadLines, iRB, unBinHash):
	Contact_bin_hash = [{},{}]#,0]
	print '\tStart binning. Contact analysis...'
	start_time = timeit.default_timer()
	for i in range(len(ReadLines)):
		parse = ReadLines[i].split()
		s0 = iRB[int(parse[0])]
		s1 = iRB[int(parse[1])]
		c = int(parse[2])

		if ( ( (s0[0] == s1[0]) and (math.fabs(s1[1] - s0[1]) > 1) ) or (s0[0] != s1[0]) ):
			if ( HashTry(unBinHash, s0) == 0 ) and ( HashTry(unBinHash, s1) == 0 ):
				#Contact_bin_hash[2] += c
				if s1[1] >= s0[1]:
					if HashTry(Contact_bin_hash[1], s0+s1) == 0: Contact_bin_hash[1][s0+s1] = c
					else: Contact_bin_hash[1][s0+s1] += c
				else:
					if HashTry(Contact_bin_hash[1], s1+s0) == 0: Contact_bin_hash[1][s1+s0] = c
					else: Contact_bin_hash[1][s1+s0] += c
				if HashTry(Contact_bin_hash[0], s0) == 0: Contact_bin_hash[0][s0] = c
				else: Contact_bin_hash[0][s0] += c
				if HashTry(Contact_bin_hash[0], s1) == 0: Contact_bin_hash[0][s1] = c
				else: Contact_bin_hash[0][s1] += c

		if i % 10000000 == 0:
			elp = timeit.default_timer() - start_time
			print '\t\tContact analyzing progress: %i, time elapsed: %.2f, memory sized: %.2f' %( i, elp, sys.getsizeof(Contact_bin_hash) )

	Keys = Contact_bin_hash[1].keys()
	for key in Keys:
		if Contact_bin_hash[1][key] == 1: 
			del Contact_bin_hash[1][key]
			Contact_bin_hash[0][key[:2]] += -1
			Contact_bin_hash[0][key[2:]] += -1
			#Contact_bin_hash[2] += -1
	
	elp = timeit.default_timer() - start_time
	print '\tContact analyzing full time: %.2f memory sized: %i' % (elp, sys.getsizeof(Contact_bin_hash) )
	return Contact_bin_hash

def IcedMatrixRead(ReadLines, iRB):
	IMR = {},{} #IMR_Hash - hash of iced matrix
	Keys = set()
	for line in ReadLines:
		parse = line.split()
		s0 = iRB[int(parse[0])]
		s1 = iRB[int(parse[1])]
		c = float(parse[2])
		IMR[0][s0+s1] = c
		IMR[0][s1+s0] = c
		if HashTry(IMR[1], s0) == 1: IMR[1][s0] += c
		else: IMR[1][s0] = c
		if HashTry(IMR[1], s1) == 1: IMR[1][s1] += c
		else: IMR[1][s1] = c
	for i in IMR[0]: IMR[0][i] = int(IMR[0][i]/IMR[1][i[:2]]*1e6)
	return IMR[0]

def iContactBinStat(Contact_bin_hash, IMR_hash, dist_max):
	l = 0
	Contact_bin_list = []
	for i in Contact_bin_hash[1]:
		if i[0] == i[2]:
			l = math.fabs(i[3] - i[1])
			if l > dist_max: l = dist_max
		else: l = dist_max
		#print l,
		if HashTry(IMR_hash, i) == 1:
			#print l
			Contact_bin_list.append (
				(i[0], i[1], i[2], i[3], Contact_bin_hash[1][i], 
				IMR_hash[i], Contact_bin_hash[0][i[:2]], Contact_bin_hash[0][i[2:]], l)
			)
			#print Contact_bin_hash[1][i], IMR_hash[i]
		else: pass

	return Contact_bin_list
	
def iDistanceHash(Contact_bin):
	Contact_VC = {}
	k = 0
	print '\tVC calculate for contacts between %i locus pair' % (len(Contact_bin))
	start_time = timeit.default_timer()
	temp = np.transpose( np.vstack( (
		Contact_bin[:,0],Contact_bin[:,1],Contact_bin[:,2],Contact_bin[:,3],Contact_bin[:,4],
		Contact_bin[:,5],
		np.sqrt( 1.0/Contact_bin[:,4] ),
		Contact_bin[:,6], Contact_bin[:,7],
		Contact_bin[:,8]
		) ) )
	keys = np.unique(temp[:,-1])
	Contact_VC = dict([(i,[]) for i in keys])
	for i in temp: Contact_VC[i[-1]].append(i)
	elp = timeit.default_timer() - start_time
	print '\tVC full time:', elp
	return Contact_VC

def iAbsoluteContacts(Contact_VC_key):
	Total_prc_list = []
	for i in range(len(Contact_VC_key)):
		disp = Contact_VC_key[i][6]
		c = Contact_VC_key[i][5]
		Total_prc_list.append( (Contact_VC_key[i][0],Contact_VC_key[i][1], Contact_VC_key[i][2], Contact_VC_key[i][3], c, int(c*(1-disp)), int(c*(1+disp)), Contact_VC_key[i][7], Contact_VC_key[i][8], Contact_VC_key[i][9]) )
	return Total_prc_list

def BinContactWriting(FilePath, prc_con, ChrInd):
	print '\tStart contact writing to', FilePath
	start_time = timeit.default_timer()
	prc_con.sort()
	elp = timeit.default_timer() - start_time
	print '\tDatabase sorting:',elp
	start_time = timeit.default_timer()
	File = open(FilePath+'.initialContacts','w')
	print >> File, 'chr1\tpos1\tchr2\tpos2\tpcr\tmin\tmax\tcov1\tcov2\tdist'
	for i in range(len(prc_con)):
		try:
			print >> File, '%s\t%i\t%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i' % (ChrInd[prc_con[i][0]],prc_con[i][1],ChrInd[prc_con[i][2]],prc_con[i][3],prc_con[i][4],prc_con[i][5],prc_con[i][6],prc_con[i][7],prc_con[i][8],prc_con[i][9])
		except IndexError:
			print 'index', prc_con[i]
		except KeyError:
			print 'key', prc_con[i]
	File.close()
	return prc_con