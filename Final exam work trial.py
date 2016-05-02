# -*- coding: utf-8 -*-
"""
Created on Mon May 02 05:53:20 2016

@author: lenovo
"""

import win32com.client

import scipy, scipy.optimize
from scipy.optimize import curve_fit
import scipy.integrate
from scipy.integrate import odeint
import math
import numpy as np
import matplotlib.pyplot as plt

xl= win32com.client.gencache.EnsureDispatch("Excel.Application")
xl.Interactive = False
wb=xl.Workbooks('Data.xlsx')
sheet=wb.Sheets('Data')

def getdata(sheet, Range):
    data= sheet.Range(Range).Value
    data=scipy.array(data)
    data=data.reshape((1,len(data)))[0]
    return data
    
x=getdata(sheet,"A2:A11")
y=getdata(sheet,"B2:B11")
z=getdata(sheet,"C2:C11")

def fit_func(x,p):
    [k1 , k2] =p
    k1=p[0]
    k2=p[1]
    return (1/(k1+k2))*(-(np.exp(-(k1+k2)*x))+k2*49.44)
  
def error(p,x,y):
    ycalc=fit_func(x,p)
    err=ycalc-y
    return err

def get_r2(x,y,ycalc):
    ymean=scipy.average(y)
    dymean2=(y-ymean)**2
    dycalc2=(y-ycalc)**2
    r2=1-sum(dycalc2)/sum(dymean2)
    return r2
    
pguess= [10000,1000]
plsq=scipy.optimize.leastsq(error, pguess, args=(x,y))
p=plsq[0]
k1=p[0]
k2=p[1]
print k1
print k2

ycalc=fit_func(x,p)
r2=get_r2(x,y,ycalc)


def concentration(z,t):
    [A,B]=z
    
    dAdt= (-k1-k2)*A+k2*49.44
    dBdt= (-k1-k2)*B+49.44*k1
    
    return [dAdt,dBdt]
    
A0=49.44
B0=0.00
z0=[A0,B0]
t=np.linspace(0,80,10)

s=odeint(concentration,z0,t)

A=s[:,0]
B=s[:,1]
soln=scipy.array([[t],[A],[B]])

print soln

plt.plot(t,A,'r')
plt.plot(t,B,'g')


    
    
    

