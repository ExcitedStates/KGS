#!/usr/bin/python
import os
import sys
from numpy import *
import pandas
import matplotlib
import matplotlib.pyplot as plt
import scipy.interpolate as ip

matplotlib.rcParams.update({'font.size': 22})
#data_dG = pandas.read_csv('energy_deltas_forward.csv')
inputFile = sys.argv[1]
data_dG = pandas.read_csv(inputFile)

outName = inputFile[:-3]
#print (data_dG)

x = data_dG['sample']
y = data_dG['dG']
#print (max(y))

#polynomial fit
coefficients = polyfit(x,y,6)
polynomial = poly1d(coefficients)
xs = arange(min(x), max(x), 0.1)
ys = polynomial(xs)

#ylim(-120,110)
# print data_dG
#rolling average
plt.plot(data_dG.dG, 'o',alpha=0.5)
ysmo = pandas.rolling_mean(data_dG.dG, window=10, min_periods=1)
plt.plot(ysmo, linewidth=3)
meanG=mean(data_dG.dG)
cumG = sum(data_dG.dG)


plt.savefig(outName+'.png',dpi=300)

#plot(xs, ys, linewidth=10)
#ylabel('dG (kcal/mol)')
#xlabel('conformational transition (reaction coordinate)')

#spline interpolation
w = ones([len(x),1])
spl = ip.splrep(x,y,w,s=10000)
xn = linspace(0,1.0,100)
sple = ip.splev(xs,spl)
#plot(xs,sple,linewidth=3)
#show()



fig = plt.figure()
ysmo_H = pandas.rolling_mean(data_dG.dH, window=10, min_periods=1)
plt.plot(ysmo_H, linewidth=3,alpha=0.5)
ysmo_S = pandas.rolling_mean(data_dG.dS, window=10, min_periods=1)
plt.plot(ysmo_S, linewidth=3,alpha=0.5)
plt.savefig('energyPlots/dU_dS_rollMean.png',dpi=300)

fig = plt.figure()
plt.plot(data_dG.dH,data_dG.dS, '.')
# plt.plot(ysmo_H,ysmo_S, '.')
meanH=mean(data_dG.dH)
meanS=mean(data_dG.dS)
cumH=sum(data_dG.dH)
cumS=sum(data_dG.dS)
# meanH=mean(ysmo_H)
# meanS=mean(ysmo_S)
# cumH=sum(ysmo_H)
# cumS=sum(ysmo_S)
plt.plot(meanH,meanS,'ro')
# plt.plot(cumH,cumS,'rs')
print "Average change: dH",meanH,"dS",meanS,"dF",meanG
print "Cumulative change: dH",cumH,"dS",cumS,"dF",cumG
plt.savefig('energyPlots/dU_vs_dS.png',dpi=300)