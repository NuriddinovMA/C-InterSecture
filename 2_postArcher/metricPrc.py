import postArcher_func as pAf
import numpy as np
import matplotlib.pyplot as plt
reload(pAf)

pair = 'HumanChicken_','HumanMouse_','MouseChicken_'
pair_out = ('hg','gh'),('hm','mh'),('mg','gm')
track_path = "//fishminja/share/ArChER/tracks/"
track_file = ("IMR90","CEF"),("IMR90","MKEF"),("MKEF","CEF")
suf = ".prc.h1.8binfr.40kb.metric.bedGraph"

spl = 230
plt.cla()
for s in range(2):
	for p in range(3):
		BEDM = pAf.readBED(track_path+pair[p]+track_file[p][s]+suf,40000,'values')
		L = []
		Keys = sorted(BEDM.keys())
		for key in Keys: L.append(BEDM[key])
		L.sort()
		ln = len(L)
		X = []
		Y = []
		k = int(L[0])
		for i in range(ln):
			if L[i] > k:
				X.append(k)
				Y.append(1.0*i/ln)
				k += 1
		spl += 1
		plt.subplot(spl)
		plt.plot(X,Y,'-')
		plt.title(pair_out[p][s], size=9)
		xl = [i for i in range(0,55,5)]
		yl = [i/10.0 for i in range(0,11)]
		plt.xlim((0,55))
		plt.xticks(xl,xl,size=6, rotation=60)
		plt.yticks(yl,yl,size=6)
		plt.grid()
plt.subplots_adjust(hspace=0.35)
plt.subplots_adjust(wspace=0.25)
plt.savefig('GO/metricstat.png', dpi = 400)
plt.clf()