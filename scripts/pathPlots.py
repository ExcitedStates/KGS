#!/usr/bin/python
import sys
import os
import numpy as np
import fnmatch
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from combineKGSPath_steps import extractPath
from math import log
import csv
from matplotlib.ticker import MultipleLocator, LinearLocator, FormatStrFormatter

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Arial']}) #Helvetica, Avant Gard, Computer Modern Sans serif
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

from pylab import rcParams
rcParams['figure.figsize'] = 6.5, 6.5*3/4

ftSize=24
mpl.rc('xtick', labelsize=ftSize) 
mpl.rc('ytick', labelsize=ftSize)

def getData(pdbPath,pathList,reversePathList):
	
	maxCons=0
	maxRevCons=0
	nullspaceDofs=0
	maxDofs=0
	revMaxDofs=0
	forwardName="init"
	targetName="target"
	
	with open("output.txt") as outputFile:
		for line in outputFile:
			if "--initial " in line:
				forwardName = line[line.rfind("/")+1:line.rfind(".pdb")]
			if "--target " in line:
				targetName = line[line.rfind("/")+1:line.rfind(".pdb")]
			if ("Total DOFs:" in line): 
				tokens = line.split(' ')
				totalDof= tokens[2]
				cycleDof = tokens[5]
				nullspaceDofs=int(tokens[9])
				maxDofs = int(totalDof[:-1]) - int(cycleDof[:-1])+nullspaceDofs
			# if ("Total remaining DOFs:" in line):
			# 	tokens = line.split(' ')
			# 	dofString = tokens[3]
			# 	maxDofs = int(dofString[:-1])
			if ("Total DOFs in target:" in line):
				tokens = line.split(' ')
				totalDof= tokens[4]
				cycleDof = tokens[7]
				if len(tokens) >= 10:
					revNullspaceDofs = int(tokens[11])
				else:
					revNullspaceDofs = nullspaceDofs
				maxDofsRev = int(totalDof[:-1]) - int(cycleDof[:-1])+revNullspaceDofs
				if maxDofsRev != maxDofs:
					print "Forward and reverse proteins have different dofs"
			# if ("Total remaining DOFs in target:" in line):
			# 	tokens = line.split(' ')
			# 	dofString = tokens[5]
			# 	maxDofsRev = int(dofString[:-1])
				# if maxDofsRev != maxDofs:
				# 	print "Forward and reverse proteins have different dofs"
	
	print "Max forward dofs: "+str(maxDofs)
	print "Max reverse dofs: "+str(maxDofsRev)
	
	pathFileSepIdx = pdbPath.rfind("/")
	outputDir = pdbPath[0:pathFileSepIdx] if pathFileSepIdx!=-1 else "."
	# sampleNum = int(pdbPath[pdbPath.rfind("_")+1:pdbPath.rfind(".pdb")])
	# modelName = str(pdbPath[pdbPath.rfind("/")+1:pdbPath.rfind("_path")])

	sample_ids = []
	rmsdIni=[]
	rmsdTarget=[]
	energy=[]
	clashFreeDofs=[]
	clashConstraints=[]
	deltaH=[]
	confId=0
	
	reverse_sample_ids = []
	reverse_rmsdIni=[]
	reverse_rmsdTarget=[]
	reverse_energy=[]
	reverse_clashFreeDofs=[]
	reverse_clashConstraints=[]
	reverse_deltaH=[]
	minDofs=999999

	# print pathList,reversePathList
	
	for sample in pathList:
		sample_id = -1
		rmsdToIni = 0.0
		rmsdToTarget = 1000.0
		clashFreeDof = 0
		clashCons = 0
		en = 0
		delta = 0
		with open(outputDir+"/"+forwardName+"_new_"+str(sample)+".pdb", "r") as pdb_file:
			for line in pdb_file:
				if "REMARK\tID" in line:
					sample_id = int(line[12:])
				if "REMARK\tRMSD to initial" in line:
					rmsdToIni = float(line[25:])
				if "REMARK\tRMSD to target" in line:
					rmsdToTarget = float(line[24:])
				# if "REMARK\tClash constraints" in line:
				# 	tokens=line.split(' ')
				# 	clashCons = int(tokens[3])
				if "REMARK\tClash free dofs" in line:
					tokens=line.split(' ')
					clashFreeDof = int(tokens[4])
					if( clashFreeDof < minDofs):
						minDofs=clashFreeDof
				if "REMARK\tVdw energy" in line:
					en = float(line[20:])
				if "REMARK\tDelta H" in line:
					delta = float(line[17:])
					break

		clashCons = maxDofs-clashFreeDof
		sample_ids.append(sample_id)
		rmsdIni.append(rmsdToIni)
		rmsdTarget.append(rmsdToTarget)
		clashConstraints.append(clashCons)
		clashFreeDofs.append(clashFreeDof)
		energy.append(en)
		deltaH.append(delta)
		if clashCons > maxCons:
			maxCons=clashCons
			confId = sample_id
		# if clashCons > 30:
		# 	print "Conf "+str(sample_id)+" has more than 30 constraints"
			
	print "Conf with highest # clashes: "+str(confId)
		
	rev=False
	sample_id = -1
	rmsdToIni = 0.0
	rmsdToTarget = 1000.0
	clashFreeDof = 0
	clashCons = 0
	en = 0
	delta=0
		
	for sample in reversePathList:
		rev=True
		sample_id = -1
		rmsdToIni = 0.0
		rmsdToTarget = 1000.0
		clashFreeDof = 0
		clashCons = 0
		en = 0
		delta=0
		with open(outputDir+"/"+targetName+"_new_"+str(sample)+".pdb", "r") as pdb_file:
			for line in pdb_file:
				if "REMARK\tID" in line:
					sample_id = int(line[12:])
				if "REMARK\tRMSD to initial" in line:
					rmsdToIni = float(line[25:])
				if "REMARK\tRMSD to target" in line:
					rmsdToTarget = float(line[24:])
				# if "REMARK\tClash constraints" in line:
				# 	tokens=line.split(' ')
				# 	clashCons = int(tokens[3])
				if "REMARK\tClash free dofs" in line:
					tokens=line.split(' ')
					clashFreeDof = int(tokens[4])
					if( clashFreeDof < minDofs):
						minDofs=clashFreeDof
				if "REMARK\tVdw energy" in line:
					en = float(line[20:])
				if "REMARK\tDelta H" in line:
					delta = float(line[17:])
					break
			
		clashCons = maxDofsRev - clashFreeDof	
		reverse_sample_ids.append(sample_id)
		reverse_rmsdIni.append(rmsdToIni)
		reverse_rmsdTarget.append(rmsdToTarget)
		reverse_clashConstraints.append(clashCons)
		reverse_clashFreeDofs.append(clashFreeDof)
		reverse_energy.append(en)
		reverse_deltaH.append(delta)
		
		if clashCons > maxRevCons:
			maxRevCons=clashCons
			revConfId = sample_id
		# if clashCons > 30:
		# 	print "Conf "+str(sample_id)+" has more than 30 constraints"
		
		if delta < -10000:
			print "Huge delta at reverse conf "+str(sample_id)
			
	if rev==True:
		print "Rev Conf with highest # clashes: "+str(revConfId)
		
	for i in range(len(clashFreeDofs)):
		if(clashConstraints[i] != maxDofs-clashFreeDofs[i]):
			print "Reduced clash cons from "+str(clashConstraints[i])+" to "+str(maxDofs-clashFreeDofs[i])
			clashConstraints[i] = maxDofs-clashFreeDofs[i]
		clashFreeDofs[i]=clashFreeDofs[i]/float(maxDofs)

	for i in range(len(reverse_clashFreeDofs)):
		if(reverse_clashConstraints[i] != maxDofs-reverse_clashFreeDofs[i]):
			#print "Reduced clash cons from "+str(reverse_clashConstraints[i])+" to "+str(maxDofs-reverse_clashFreeDofs[i])
			reverse_clashConstraints[i] = maxDofs-reverse_clashFreeDofs[i]
		reverse_clashFreeDofs[i]=reverse_clashFreeDofs[i]/float(maxDofs)
		
	samples=np.arange(0,len(energy))
	# plotEnthalpyAndEntropy(samples,energy,clashConstraints,"vdW energy","independent clash constraints")
	
	return energy, clashConstraints, deltaH, clashFreeDofs, reverse_energy, reverse_clashConstraints, reverse_deltaH, reverse_clashFreeDofs
	
def plotEnthalpyAndEntropy(xData,samples1, samples2,yString1,yString2,saveFileName=None,switch=None):

	import matplotlib.pyplot as plt
	import numpy as py
	import math
	from matplotlib.ticker import MultipleLocator, LinearLocator, FormatStrFormatter

	from matplotlib import rc
	rc('font',**{'family':'sans-serif','sans-serif':['Arial']}) #Helvetica, Avant Gard, Computer Modern Sans serif
	## for Palatino and other serif fonts use:
	#rc('font',**{'family':'serif','serif':['Palatino']})
	rc('text', usetex=True)
	import math

	from pylab import rcParams
	# rcParams['figure.figsize'] = 6.5, 6.5*3/4

	ftSize=24
	mpl.rc('xtick', labelsize=ftSize) 
	mpl.rc('ytick', labelsize=ftSize)

	sizeFactor1=500
	sizeFactor2=500
	
	maxLength = len(xData)
	fig=plt.figure()
	ax = fig.add_axes([0.2,0.2,0.68,0.75])
	ax.set_xlabel('transition path [100 steps]',fontsize=ftSize)
	majorTicks = np.arange(0,maxLength,500)
	majorTickLabels = np.arange(0,(maxLength)/100+1,5)

	# minY = math.floor(min(samples1)/sizeFactor1 )*sizeFactor1/100-sizeFactor1/100
	# maxY = math.ceil(max(samples1)/sizeFactor1 )*sizeFactor1/100+sizeFactor1/100
	minY=math.floor(min(samples1))
	maxY = math.ceil(max(samples1))
	# stepWidth = math.ceil((maxY - minY)/5)

	# majorYTickLabels = np.arange(minY,maxY+1,sizeFactor1/100)
	majorLocator = MultipleLocator(sizeFactor1)
	majorFormatter = FormatStrFormatter('%d')
	# # majorYTickLabels = np.arange(minY,maxY+1,5)
	# majorLocator = MultipleLocator(stepWidth)
	
	# maxY2 = math.ceil(max(samples2))
	# stepWidth2 = math.ceil((maxY2)/4)
	majorYLocator = MultipleLocator(2)
	majorYFormatter = FormatStrFormatter('%d')
	
	# majorY2Locator = MultipleLocator(sizeFactor2/100)
	majorY2Locator = MultipleLocator(5)
	majorY2Formatter = FormatStrFormatter('%d')
	# minorY2Locator = MultipleLocator(1)
	
	# for the minor ticks, use no labels; default NullFormatter
	# ax.xaxis.set_minor_locator(minorLocator)
	ax.set_xticks(majorTicks)
	ax.set_xticklabels(majorTickLabels)
	# ax.yaxis.set_major_locator(majorYLocator)
	ax.yaxis.set_major_formatter(majorYFormatter)
	# ax.yaxis.set_minor_locator(minorLocator)
	
	# ax.set_yticklabels(majorYTickLabels)
	# ax.set_yticklabels(ax.get_yticklabels())
	# for the minor ticks, use no labels; default NullFormatter
	# ax.yaxis.set_minor_locator(minorLocator)

	ax.plot(xData,samples1, 'b-')

	# yminLine,ymaxLine = ax.get_ylim();
	# if switch:
	# 	ax.plot([xData[switch], xData[switch]], [yminLine,max(samples1[switch],samples1[switch+1])],'k--')  

	# if switch: 
	# 	ax.plot([xData[switch], xData[switch]], [mayorYTickLabels[0],samples1[switch]],'k--')  
	
	# Make the y-axis label and tick labels match the line color.
	ax.set_ylabel(yString1, color='b',fontsize=ftSize)
	for tl in ax.get_yticklabels():
		tl.set_color('b')
	
	ax2 = ax.twinx()
	ax2.plot(xData,samples2, 'r.')
	# ax2.yaxis.set_major_locator(majorY2Locator)
	ax2.yaxis.set_major_formatter(majorY2Formatter)
	# for the minor ticks, use no labels; default NullFormatter
	# ax2.yaxis.set_minor_locator(minorY2Locator)

	ax2.set_ylabel(yString2, color='r',fontsize=ftSize)
	for tl in ax2.get_yticklabels():
		tl.set_color('r')
		
	yminLine,ymaxLine = ax2.get_ylim();
	if switch:
		ax2.plot([xData[switch], xData[switch]], [yminLine,ymaxLine],'k--')    
		# ax2.plot([xData[switch], xData[switch]], [0,12],'k--')                                      
		
	# if saveFileName == "vdwAndClashes_completePath.png":
	# 	ax.text(50,-7700,'remaining energy barrier',transform=ax.transData,fontsize=ftSize)
	# 	ax.annotate('', xy=(2000,-7900),xytext=(1200,-7800),arrowprops=dict(arrowstyle="->",facecolor='black'),transform=ax.transData)  

	if saveFileName:
		plt.savefig(saveFileName,transparent='True',format='png',dpi=308)
	plt.clf()
	plt.close()
	
def plotEnthalpy(xData,samples1,yString1,saveFileName=None,switch=None):

	import matplotlib.pyplot as plt
	import numpy as py
	import math
	from matplotlib.ticker import MultipleLocator, LinearLocator, FormatStrFormatter
	from matplotlib import rc
	rc('font',**{'family':'sans-serif','sans-serif':['Arial']}) #Helvetica, Avant Gard, Computer Modern Sans serif
	## for Palatino and other serif fonts use:
	#rc('font',**{'family':'serif','serif':['Palatino']})
	rc('text', usetex=True)

	from pylab import rcParams
	rcParams['figure.figsize'] = 6.5, 6.5*3/4

	ftSize=24
	mpl.rc('xtick', labelsize=ftSize) 
	mpl.rc('ytick', labelsize=ftSize)

	sizeFactor1=500
	sizeFactor2=500
	
	maxLength = len(xData)
	fig=plt.figure()
	ax = fig.add_axes([0.2,0.2,0.68,0.75])
	ax.set_xlabel('transition path [100 steps]',fontsize=ftSize)
	majorTicks = np.arange(0,maxLength,500)
	majorTickLabels = np.arange(0,(maxLength)/100+1,5)
	
	ax.set_xticks(majorTicks)
	ax.set_xticklabels(majorTickLabels)

	ax.plot(xData,samples1)

	if saveFileName:
		plt.savefig(saveFileName,transparent='True',format='png',dpi=600,figsize=(6.5, 6.5*3/4))
	plt.clf()
	plt.close()


def main():
	"""
	Checks arguments and writes free energy along the tree-path of a kgs pdb-file
	"""

	if len(sys.argv)<2:
		print "Usage: "+sys.argv[0]+" <path.pdb files in a row> "
		sys.exit(1)

	samples=[]
	energy=[]
	clashConstraints=[]
	deltaH=[]
	relDofs=[]
	rev_energy=[]
	rev_clashConstraints=[]
	rev_deltaH=[]
	rev_relDofs=[]
	
	# pdbPath = sys.argv[1]
	# numRounds= int(sys.argv[2])
	
	for i in range(len(sys.argv)-1):
		pdbPath=sys.argv[i+1]
		pathList, reversePathList = extractPath(pdbPath)
		# print len(pathList),len(reversePathList)
		en, clashCons, delta, relDof, rev_en, rev_clashCons, rev_delta, rev_relDof = getData(pdbPath,pathList, reversePathList)
		energy.extend(en)
		clashConstraints.extend(clashCons)
		deltaH.extend(delta)
		relDofs.extend(relDof)
		rev_energy.extend(rev_en)
		rev_clashConstraints.extend(rev_clashCons)
		rev_deltaH.extend(rev_delta)
		rev_relDofs.extend(rev_relDof)
	# for i in range(numRounds-1):
	# 	pathList, reversePathList = extractPath(pdbPath)
	# 	en, clashCons = getData(pdbPath,pathList, reversePathList)
	# 	energy.extend(en)
	# 	clashConstraints.extend(clashCons)
	# 	print os.getcwd()
	# 	os.chdir("*round"+str(i+2))
	# 	os.chdir("output")
	# 	pdbPath="*path.pdb"
	# 	
	# os.chdir("../..")
		
	samples = np.arange(0,len(energy))
	
	deltaS=[log(y) for y in relDofs]
	Tfac=1000
	deltaG=[x-Tfac*y for x,y in zip(deltaH,deltaS)]
	
	rev_deltaS=[log(y) for y in rev_relDofs]
	rev_Tfac=Tfac
	rev_deltaG=[x-Tfac*y for x,y in zip(rev_deltaH,rev_deltaS)]

	
	d="energyPlots"
	if not os.path.exists(d):
		os.makedirs(d)
	os.chdir(d)   
	
	orig_stdout = sys.stdout
	fout = file('energy_deltas_forward.txt', 'w')
	sys.stdout = fout

	for i in range(len(samples)):
		print samples[i], deltaG[i], deltaH[i], deltaS[i]
		
	sys.stdout = orig_stdout
	fout.close()
		
	Rev=len(rev_energy)>0
	
	if Rev:
		samples = np.arange(0,len(rev_energy))			
		orig_stdout = sys.stdout
		fout = file('energy_deltas_reverse.txt', 'w')
		sys.stdout = fout

		for i in range(len(samples)):
			print samples[i], rev_deltaG[i], rev_deltaH[i], rev_deltaS[i]
			
		sys.stdout = orig_stdout
		fout.close()
		
		#Combined path
		samples = np.arange(0,len(energy)+len(rev_energy))
		rev_energy.reverse()
		rev_clashConstraints.reverse()
		rev_deltaH.reverse()
		rev_deltaS.reverse()
		rev_deltaG.reverse()
		
		switch=len(deltaG)
		
		deltaG.extend(rev_deltaG)
		energy.extend(rev_energy)
		clashConstraints.extend(rev_clashConstraints)
		deltaH.extend(rev_deltaH)
		deltaS.extend(rev_deltaS)
		
		plotEnthalpyAndEntropy(samples,energy,clashConstraints,"vdW energy [kcal/mol]","independent clash constraints","vdwAndClashes_completePath.png",switch)
		plotEnthalpy(samples,energy,"vdW energy [kcal/mol]","vdw_completePath.png")		
	else:
		plotEnthalpyAndEntropy(samples,energy,clashConstraints,"vdW energy [kcal/mol]","independent clash constraints","vdwAndClashes_forwardPath.png")
		plotEnthalpy(samples,energy,"vdW energy [kcal/mol]","vdw_forwardPath.png")		
	orig_stdout = sys.stdout
	# fout = file('energy_deltas_complete.txt', 'w')
	# sys.stdout = fout
	# 
	# for i in range(len(samples)):
	# 	print samples[i], deltaG[i], deltaH[i], deltaS[i]
	# 	
	# sys.stdout = orig_stdout
	# fout.close()
	with open('energy_deltas_complete.csv','wb') as csvfile:
		spamwriter = csv.writer(csvfile, delimiter=',')
		spamwriter.writerow(['sample','dG','dH','dS'])
		for i in range(len(samples)):
			spamwriter.writerow([samples[i], deltaG[i], deltaH[i], deltaS[i]])
	
if __name__ == "__main__":
	main()
