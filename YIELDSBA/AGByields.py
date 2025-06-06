#! /usr/bin/python
import pylab as pl
import matplotlib.pyplot as plt
import pickle
import math
from scipy import*
import array
import string,re,sys,os
import numpy as np

data = np.loadtxt('yieldsbaGallo.dat',skiprows=1)

plt.plot(data[:,0],data[:,1]*1.5,'k--',lw=4)
plt.plot(data[:,0],data[:,2]*3.,'r--',lw=4)

data = np.loadtxt('Cristallo1.3.dat',skiprows=1)

plt.plot(data[:,1],data[:,2]+data[:,3]+data[:,4]+data[:,5]+data[:,6]+data[:,7]+data[:,8],'b-',lw=2)

data = np.loadtxt('Cristallo1.5.dat',skiprows=1)

plt.plot(data[:,1],data[:,2]+data[:,3]+data[:,4]+data[:,5]+data[:,6]+data[:,7]+data[:,8],'r-')


data = np.loadtxt('Cristallo2.0.dat',skiprows=1)

plt.plot(data[:,1],data[:,2]+data[:,3]+data[:,4]+data[:,5]+data[:,6]+data[:,7]+data[:,8],'g-')

data = np.loadtxt('Cristallo2.5.dat',skiprows=1)

plt.plot(data[:,1],data[:,2]+data[:,3]+data[:,4]+data[:,5]+data[:,6]+data[:,7]+data[:,8],'m-')

data = np.loadtxt('Cristallo3.0.dat',skiprows=1)

plt.plot(data[:,1],data[:,2]+data[:,3]+data[:,4]+data[:,5]+data[:,6]+data[:,7]+data[:,8],'k-')



plt.axis([0,0.04,1.e-10,1.e-6])

plt.show()
