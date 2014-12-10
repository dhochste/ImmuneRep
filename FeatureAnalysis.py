import re
import matplotlib.pyplot as plt
import numpy as np

# Frequency Analysis
def FrequencyAnalysis(samplename,dirname='../data/',segment='V',plotflag=0):
	freqs=dict()
	if segment=='V':
		filename=dirname+samplename+'.V.germdata.txt'
		seglist=[re.findall(r'IGHV[0-9]*-[0-9]+/*[0-9]*',line) for line in open(filename)]
		for i in seglist:
			for j in i:
				if j in freqs:
					freqs[j]+=1
				else:
					freqs[j]=1
	nsamp=sum([freqs[k] for k in freqs.keys()])
	for k in freqs.keys():
		freqs[k]=float(freqs[k])/float(nsamp)
	# print [list(re.findall(r'IGHV([0-9]*)-([0-9]+)',k)[0]) for k in freqs.keys()]
	# segnums=[i for i in re.findall(r'IGHV([0-9]*)-([0-9]+)',k) for k in freqs.keys()]
	segnames=freqs.keys()
	segnums=[[int(i) for i in list(re.findall(r'IGHV([0-9]*)-([0-9]+)',k)[0])] for k in freqs.keys()]
	segnames=[i for j,i in sorted(zip(segnums,segnames))]
	print segnames

	if plotflag:
		ind=np.arange(len(freqs.keys()))
		fs=[freqs[k] for k in segnames]
		plt.bar(ind,fs)
		plt.xticks(ind+.4,segnames,rotation='vertical')
		plt.show()
	return freqs

def FreqDictCombiner(freqsdicts):
	segnames=[]
	for fdict in freqsdicts:
		segnames=list(set(segnames) | set(fdict.keys()))



samplename='SRR1298742'
# segnames,frequencies=FrequencyAnalysis(samplename)
freqs=FrequencyAnalysis(samplename,plotflag=1)
# print freqs.keys()
# print [freqs[k] for k in freqs.keys()]
# plt.bar(range(len(freqs.keys())),[freqs[k] for k in freqs.keys()])
