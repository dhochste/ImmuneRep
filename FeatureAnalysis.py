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

def FrequencyCorrelation(freqarray,segnames,plotflag=0):
	corr=np.corrcoef(freqarray)
	if plotflag==1:
		plt.pcolor(corr,cmap='RdBu',vmin=0,vmax=1)
		plt.xticks(np.arange(float(len(corr)))+.5,samplenames)
		plt.yticks(np.arange(float(len(corr)))+.5,samplenames,rotation='vertical')
		plt.show()
	return corr


# this function doesn't work or do anything yet
def FrequencyCluster(samplenames,dirname='../data/',segment='V'):
	freqs=FrequencyCountN(samplenames,dirname=dirname,segment=segment)
	freqarray,segnames=FrequencyArray(freqs)
	model=AgglomerativeClustering(linkage='ward')
	model.fit(freqarray)
	print model.labels_


class FrequencyAnalyzer:
	def __init__(self,samplenames,dirname='../data/',segment='V'):
		self.samplenames=samplenames
		self.dirname=dirname
		self.segment=segment

	def count(self,plotflag=0):
		self.freqs=FrequencyCountN(samplenames=self.samplenames,dirname=self.dirname,segment=self.segment,plotflag=plotflag)
		return self.freqs

	def sort(self):
		self.segnames=SortNames(self.freqs)
		return self.segnames

	def asarray(self):
		if not hasattr(self,'freqs'):
			self.count()
		self.freqarray,self.segnames=FrequencyArray(freqs=self.freqs)
		return self.freqarray

	def corr(self,plotflag=0):
		if not hasattr(self,'freqarray'):
			self.asarray()
		self.corrcoef=FrequencyCorrelation(freqarray=self.freqarray,segnames=self.segnames,plotflag=plotflag)
		return self.corrcoef



############ Clustering #################
from Bio.Seq import Seq
from Bio import motifs
from Bio.Alphabet import IUPAC
import difflib


samplename='SRR1298742'
dirname='../data/'
plotflag=0

filename=dirname+samplename+'.VDJ.H3.L3.CH1.fa'
reseqname=samplename+r'\.([0-9]*)'
reV=r'(IGHV[0-9]*-[0-9]+/*[0-9]*-*[0-9]*)'
reJ=r'(IGHJ[0-9]*)'
reseqs=r'([A-Z]*)'
reall=reseqname+r'.*'+reV+r'.*'+reJ+r'\s[0-9]*\s[0-9]*;'+reseqs
# reall=reseqname+r'.*'+reV+r'.*'+reJ+r'\s[0-9]*\s[0-9]*;'
# print re.findall(reall,'>SRR1298742.1 IBP32IW01DTCHY length=506;IGHV1-2 292 17;;IGHJ6 46 1;RAKGASDSNYAGGMDVW;;IGHG 27 0;;YTFSGYYMH;GWINPNSGGTNYA;;')
seqnames,Vs,Js,seqs=zip(*[re.findall(reall,line)[0] for line in open(filename) if line[0]=='>'])
seqVJ=dict()
for seqname,V,J,seq in zip(seqnames,Vs,Js,seqs):
	if (V,J) in seqVJ:
		seqVJ[(V,J)]['seqnames'].append(seqname)
		seqVJ[(V,J)]['seqs'].append(seq)
		seqVJ[(V,J)]['count']+=1
	else:
		seqVJ[(V,J)]={'seqnames':[seqname],'seqs':[seq],'count':1}

VJ=seqVJ.keys()
segnum=[re.findall(r'IGHV([0-9]*)-([0-9]+)/*[0-9]*-*[0-9]*.*IGHJ([0-9]*)',v+j)[0] for v,j in VJ]
VJ=[j for i,j in sorted(zip(segnum,VJ))]
counts=[seqVJ[vj]['count'] for vj in VJ]

if plotflag==1:
	ind=np.arange(len(VJ))
	plt.bar(ind,counts)
	plt.xticks(ind+.4,[v+' '+j for v,j in VJ],rotation='vertical')
	ax=plt.gca()
	for tick in ax.xaxis.get_major_ticks():
		tick.label.set_fontsize(8) 
	

	plt.figure()
	plt.hist(counts,50)
	plt.xlabel('occurences of VJ pair')
	plt.show()

def PSSM(seqs):
	n=max([len(seq) for seq in seqs])
	seqalign=[Seq(seq+'X'*(n-len(seq)),alphabet=IUPAC.ExtendedIUPACProtein) for seq in seqs]
	M=motifs.create(seqalign)
	sm=difflib.SequenceMatcher()
	ratios=[]
	for seq in seqalign:
		sm.set_seqs(M.consensus,seq)
		ratios.append((1.-sm.ratio())*float(len(M.consensus)))
	print M
	plt.hist(ratios,20)
	plt.show()

def CountClones(seqs):
	n=max([len(seq) for seq in seqs])
	seqalign=[seq+'X'*(n-len(seq)) for seq in seqs]
	clones=list(set(seqalign))
	counts=[]
	for clone in clones:
		counts.append(seqalign.count(clone))

	ind=np.arange(len(counts))
	plt.bar(ind,counts)
	plt.xticks(ind+.4,clones,rotation='vertical')
	plt.show()







# PSSM(seqVJ[('IGHV1-2','IGHJ6')]['seqs'])
# PSSM(seqVJ[('IGHV1-8','IGHJ5')]['seqs'])
# PSSM(seqVJ[('IGHV1-69','IGHJ6')]['seqs'])
# PSSM(seqVJ[('IGHV4-34','IGHJ6')]['seqs'])

CountClones(seqVJ[('IGHV4-34','IGHJ6')]['seqs'])


################################################
# samplenames=['SRR1298742','SRR1298742']
# segnames,frequencies=FrequencyCount(samplename)
# freqs=FrequencyCountN(samplenames,segment='V',plotflag=1)
# print corr
