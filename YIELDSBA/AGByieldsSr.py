#! /usr/bin/python
import pylab as pl
import matplotlib.pyplot as plt
import pickle
import math
from scipy import*
import array
import string,re,sys,os
import numpy as np

data = np.loadtxt('CristalloSr.dat',skiprows=1)
data2 = np.loadtxt('CristalloBa.dat',skiprows=1)
Sr=[]
Ba=[]
#for i in range (0,4):
#    for j in range (0,7):

for i in range (0,4):
    Sr=data[i*7:(i+1)*7,2]+data[i*7:(i+1)*7,3]+data[i*7:(i+1)*7,4]+data[i*7:(i+1)*7,5]+data[i*7:(i+1)*7,6]+data[i*7:(i+1)*7,7]+data[i*7:(i+1)*7,8]
    Ba=data2[i*7:(i+1)*7,2]+data2[i*7:(i+1)*7,3]+data2[i*7:(i+1)*7,4]+data2[i*7:(i+1)*7,5]+data2[i*7:(i+1)*7,6]+data2[i*7:(i+1)*7,7]+data2[i*7:(i+1)*7,8]
    
    SrBa=log10(Sr/Ba)  -4.889 + 4.334

    plt.plot(log10(data[i*7:(i+1)*7,1]/0.012),SrBa,lw=2)


#
#plt.plot(data[7:14,1],data[7:14,2]+data[7:14,3]+data[7:14,4]+data[7:14,5]+data[7:14,6]+data[7:14,7]+data[7:14,8],'r-')
#
#plt.plot(data[14:21,1],data[14:21,2]+data[14:21,3]+data[14:21,4]+data[14:21,5]+data[14:21,6]+data[14:21,7]+data[14:21,8],'g-')
#
#plt.plot(data[21:28,1],data[21:28,2]+data[21:28,3]+data[21:28,4]+data[21:28,5]+data[21:28,6]+data[21:28,7]+data[21:28,8],'m-')
#
#
#plt.plot(data[28:35,1],data[28:35,2]+data[28:35,3]+data[28:35,4]+data[28:35,5]+data[28:35,6]+data[28:35,7]+data[28:35,8],'k-')
#

plt.xlabel('[M/H]')
plt.ylabel('[Sr/Ba]')
plt.axis([-2,0.4,-1.,1.])

plt.show()
