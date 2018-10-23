import postArcher_func as pAf
import numpy as np
import matplotlib.pyplot as plt

reload(pAf)

path_track = "//fishminja/share/ArChER/tracks/old/"
file_track = "HumanMouse_IMR90.prc.h1.8binfr.40kb.metric.bedGraph"
path_features = "C:/Desktop/Arbeit/AllData/Genome/CNE/"
file_features = "HAR_hg38.bed"
f = 50
V = pAf.trackValues(pAf.readBED(path_track+file_track,40000,'values'))
F = pAf.trackFeatures(pAf.readBED(path_features+file_features,40000,'values'))
c = 1.0/25
n = 0
#for key in F:
	#FT = {key:F[key]}
	#D = pAf.distributionCalc(FT,V,f,100)
	#plt.fill_between(range(2*f+1),D[1]+3*D[2],D[1]-3*D[2],color='0.8')
	#plt.plot(D[1],color=(0,1 - n,0 + n))
	#plt.plot(D[0],color=(1,0 + n,0))
	#n += c
D = pAf.distributionCalc(F,V,f,100)
plt.fill_between(range(2*f+1),D[1]+3*D[2],D[1]-3*D[2],color='0.9')
plt.plot(D[1],color='0.2')
plt.plot(D[0],color='b')
plt.xticks(range(0,2*f+1,10),range(-1*f,f+1,10))
plt.show()
