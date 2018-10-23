import postArcher_func as pAf
import numpy as np
reload(pAf)

pair = 'HumanChicken_','HumanMouse_','MouseChicken_'
pair_out = 'hg','hm','mg'
gene_path = "C:/Desktop/Arbeit/AllData/Genome/RefGene/"
track_path = "//fishminja/share/ArChER/tracks/"
gene_file = "hg38","hg38","mm10"
suf = '.refGene.bed'
track_file = "IMR90.prc.h1.8binfr.40kb.metric.bedGraph","IMR90.prc.h1.8binfr.40kb.metric.bedGraph","MKEF.prc.h1.8binfr.40kb.metric.bedGraph"

for s in range(3):
	BEDM = pAf.readBED(track_path+pair[s]+track_file[s],40000,'values')
	BEDG = pAf.readBED(gene_path+gene_file[s]+suf,40000,'features')

	R = set([])
	Keys = sorted(BEDG.keys())
	for key in Keys:
		for i in BEDG[key]: R.add(i)
	R = list(R)
	f = open('GO/genes.%s.all.list' % gene_file[s],'w')
	for j in R: print >> f, j
	f.close()

	L = []
	Keys = sorted(BEDM.keys())
	for key in Keys: L.append(BEDM[key])
	L.sort()
	ln = len(L)
	n = [33,25,20,]
	P = []
	for i in n: P.append( L[int((100-i)*ln/100.0)] )
	print P


	for i in range(len(P)):
		L = set([])
		BEDF = pAf.filterBED(BEDM,BEDG,P[i],8)
		Keys = sorted(BEDF.keys())
		f = open('GO/%s/genes.%s.%i.bed' % (pair_out[s],pair_out[s],n[i]),'w')
		for key in Keys:
			for j in BEDF[key]: 
				L.add(j)
				print >> f, '%s\t%i\t%i\t%s' % (key[0],key[1]*40000,key[1]*40000 + 39999, j)
		f.close()
		f = open('GO/%s/genes.%s.%i.list' % (pair_out[s],pair_out[s],n[i]),'w')
		for j in L: print >> f, j
		f.close()
		f = open('GO/%s/genes.%s.%i.random.list' % (pair_out[s],pair_out[s],n[i]),'w')
		L = np.random.choice(R,len(L),replace=False)
		for j in L: print >> f, j
		f.close()