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
	M = {}
	for chr1 in parsedNet:
		M[chr1] = {}
		for i in parsedNet[chr1]:
			for j in i:
				start1 = j[1]
				end1 = start1+j[2]
				chr2 = j[3]
				dir = j[6]
				if j[0] == 'fill':
					start2 = j[4]
					end2 = j[4] + j[5] - 1
					if dir == 1: pass
					else: start2,end2 = end2,start2
					key = (chr1,start1,end1,chr2,start2,end2,dir)
					M[chr1][key] = [(start1,start1,start2,start2),]
				else:
					start2 = j[4] - 1*dir
					end2 = j[4] + j[5]*dir
					M[chr1][key].append((start1,end1,start2,end2))
	for chrName in M:
		for key in M[chrName]:
			m = M[chrName][key]
			temp = []
			if len(m) == 1: M[chrName][key] = [(key[1],key[2],key[4],key[5])]
			else:
				for i in range(1,len(m)): temp.append( (m[i-1][1],m[i][0],m[i-1][3],m[i][2]) )
				temp.append( (m[i][1],key[2],m[i][3],key[5]) )
				M[chrName][key] = temp
	color = colorList()
	f1 = open(name+'.pre.mark', 'w')
	f2 = open(name+'.2D.ann', 'w')
	print >> f1, 'chr1\tstart1\tend1\tchr2\tstart2\tend2'
	print >> f2, 'chr1\tstart1\tend1\tchr2\tstart2\tend2\tcolor\tcomment'
	chrNames = sorted(M.keys())
	for chrName in chrNames:
		syn = sorted(M[chrName].keys())
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
				for i in M[chrName][key]: print >> f1, '%s\t%i\t%i\t%s\t%i\t%i' % (key[0],i[0],i[1],key[3],i[2],i[3])
			else: pass
	f1.close()
	f2.close()
	return M


def pre2mark(M,name):
	M2 = []
	chrNames = sorted(M.keys())
	
	for chrName in chrNames:
		syn = sorted(M[chrName].keys())
		for s in syn:
			ln = len(M[chrName][s])
			coor = M[chrName][s][0]
			mark = [s[0],coor[0],coor[1],s[3],coor[2],coor[3]]
			gap = 0,0
			for i in range(1,ln):
				coor = M[chrName][s][i]
				gap = coor[0]-mark[2],coor[2]-mark[5]
				if gap[0] < 200 and gap[1] < 200: mark[2],mark[5] = coor[1],coor[3]
				else: pass
				if gap[0] > 200 or gap[1] > 200 or i == ln-1:
					c = min(mark[5]-mark[4],mark[2]-mark[1])/150
					if c > 1:
						c1,c2 = (mark[2]-mark[1])/c, (mark[5]-mark[4])/c
						for k in range(c-1): M2.append( [s[0],mark[1]+k*c1,mark[1]+k*c1+c1,s[3],mark[4]+k*c2,mark[4]+k*c2+c2] )
						k += 1
						M2.append( [s[0],mark[1]+k*c1,mark[2],s[3],mark[4]+k*c2,mark[5]] )
					else: M2.append( [s[0],mark[1],mark[2],s[3],mark[4],mark[5]] )
					mark = [ s[0],coor[0],coor[1],s[3],coor[2],coor[3] ]
				else: pass
	f1 = open(name+'.mark', 'w')
	print >> f1, 'chr1\tstart1\tend1\tchr2\tstart2\tend2'
	for m in M2:
		if len(m[0]) < 6 and len(m[3]) < 6: print >> f1, '%s\t%i\t%i\t%s\t%i\t%i' % (m[0],m[1],m[2],m[3],m[4],m[5])
		else: pass
	f1.close()
	del M2