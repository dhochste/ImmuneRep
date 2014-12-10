import re
import matplotlib.pyplot as plt
import numpy as np
from sklearn.cluster import AgglomerativeClustering



############ Frequency Analysis ###################

# utility function that pulls out the segment names from the dict and sorts them
def SortNames(freq):
	segnames=freq.keys()
	segment=re.findall(r'IGH(.)',segnames[0])[0]
	if segment in ['V','D']:
		numre=r'IGH{0}([0-9]*)-([0-9]+)/*[0-9]*-*[0-9]*'.format(segment)
	elif segment=='J':
		numre=r'IGHJ([0-9]*)'

	segnums=[[int(i) for i in list(re.findall(numre,k)[0])] for k in freq.keys()]
	segnames=[i for j,i in sorted(zip(segnums,segnames))]
	return segnames

#utility function that transforms the dict of frequencies into an array
def FrequencyArray(freqs):
	segnames=SortNames(freqs)
	freqarray=np.array([freqs[k] for k in segnames]).T
	return freqarray,segnames


#main processing function that counts segment name occurences
def FrequencyCount(samplename,dirname='../data/',segment='V',plotflag=0):
	# create regular expressions for parsing the file based on the type of segment
	if segment in ['V','D']:
		namere=r'IGH{0}[0-9]*-[0-9]+/*[0-9]*-*[0-9]*'.format(segment)
		numre=r'IGH{0}([0-9]*)-([0-9]+)/*[0-9]*-*[0-9]*'.format(segment)
	elif segment=='J':
		namere=r'IGHJ[0-9]*'
		numre=r'IGHJ([0-9]*)'

	# look through germdata to find segment names
	filename=dirname+samplename+'.'+segment+'.germdata.txt'
	seglist=[re.findall(namere,line) for line in open(filename)]
	
	# create a dictionary of counts for each segment name and normalize to get frequency
	freq=dict()
	for i in seglist:
		for j in i:
			if j in freq:
				freq[j]+=1
			else:
				freq[j]=1
	nsamp=sum([freq[k] for k in freq.keys()])
	for k in freq.keys():
		freq[k]=float(freq[k])/float(nsamp)

	segnames=SortNames(freq)

	# plot the frequencies on a bar plot
	if plotflag:
		ind=np.arange(len(freq.keys()))
		fs=[freq[k] for k in segnames]
		plt.bar(ind,fs)
		plt.xticks(ind+.4,segnames,rotation='vertical')
		plt.show()
	return freq

# aggregating function that runs the processing on multiple samples
def FrequencyCountN(samplenames,dirname='../data/',segment='V',plotflag=0):
	freqlist=[FrequencyCount(samplename,dirname=dirname,segment=segment) for samplename in samplenames]
	
	# find the union of all the lists of segment names
	segnames=[]
	for freq in freqlist:
		segnames=list(set(segnames) | set(freq.keys()))

	# list frequencies for each sample for each segment name (zero if sample not observed)
	freqs=dict()
	for k in segnames:
		freqs[k]=[]
		for freq in freqlist:
			if k in freq:
				freqs[k].append(freq[k])
			else:
				freqs[k].append(0.)

	if plotflag==1:
		freqarray,segnames=FrequencyArray(freqs)
		ind=np.arange(len(segnames))
		colors=plt.cm.spectral(np.linspace(0,1,freqarray.shape[0]))

		barwidth=0.8/freqarray.shape[0]

		for i,(fs,color) in enumerate(zip(freqarray,colors)):
			plt.bar(ind+barwidth*i,fs,barwidth,color=color)
			
		plt.xticks(ind+barwidth,segnames,rotation='vertical')
		plt.show()

	return freqs

# this function doesn't work or do anything yet
def FrequencyCluster(samplenames,dirname='../data/',segment='V'):
	freqs=FrequencyCountN(samplenames,dirname=dirname,segment=segment)
	freqarray,segnames=FrequencyArray(freqs)
	model=AgglomerativeClustering(linkage='ward')
	model.fit(freqarray)
	print model.labels_





samplenames=['SRR1298742','SRR1298742']
# segnames,frequencies=FrequencyCount(samplename)
# freqs=FrequencyCountN(samplenames,segment='V',plotflag=1)
FrequencyCluster(samplenames,segment='V')
