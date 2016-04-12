import scipy
import numpy as npy
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import scipy.integrate

Pma=500*10**-10 # permeability of component A in cm^3(STP)/scm^2cmHg
a=10.0
Pmb=10*Pma # permeability of component B in cm^3(STP)/scm^2cmHg
nf=1*10**6 # number of moles in feed in cm^3(STP)/s
theta= 0.2
# theta= np/nf= no. of moles in permeate/ no. of moles in feed
np=theta*nf
nr=nf-np
Pp=19.0 # permeate pressure in cmHg
Pf=190.0 # feed pressure in cmHg
r=Pp/Pf # pressure ratio
xfa=0.21 # mole fraction of oxygen in the feed
L=2.54*10**-3 # lenght of the membrane

np=theta*nf
nr=nf-np

# we assume completely mixed model to get the initial guess values to solvethe differental equations


def equation(yp):
    f=yp/(yp-1)-(a*(xfa-r*yp))/((1-xfa)-r*(1-yp))
    return f
A=fsolve(equation,0.3)

print A # initial permeate composition at x=xfa
yp= A[0]
print yp


npa=theta*nf*yp # rate of transger of component A
Na=(Pma/L)*(xfa*Pf-yp*Pp) # flux of transfer of A
Am=npa/Na # area of the membrane
print Am*10**-4 # area required in m^2
Pmb1=Pmb/L
Kf=(nf/Pmb1*Am*Pp)*(10**-6*10**-4*10**-4/10**-6*10**-4*10**-2)


k=input("0-cocurrent 1-countercurrent")

if (k==0):
    
    " for cocurrent flow "
    
    
    def deriv(z,A):
        [x,y]=z
        
        dxdA=(1/Kf)*((x-y)/(y-xfa))*(a*(1-x)*(x*r-y)-x*((1-x)*r)-(1-y))
        dydA=(1/Kf)*((x-y)/(x-xfa))*(a*(1-y)*(x*r-y)-y*((1-x)*r)-(1-y))
        
        # reference: A simple analysis for gas separation membranes by Richard A. Davis and Orville E. Sandall
        # university of Minnesota Duluth
        return [dxdA,dydA]
    
    # initial conditions
    x0=xfa+0.00001
    y0=yp
    z0=[x0,y0]
    print z0
    # time grid for integration
    A=scipy.linspace(0.0,1.0,11.0)
    print A
    
    p=scipy.integrate.odeint(deriv, z0, A)
    
    x=p[:,0]
    y=p[:,1]
    soln=scipy.array([[A],[x],[y]])
    
    print soln
    
    plt.plot(A,x,'g')
    plt.show()
    plt.plot(A,y,'b')
    plt.show()

if (k==1):
    "for countercurrent flow "
    
    
    def f(x,y):
        f1=y/(y-1)-(a*(x-r*y))/((1-x)-r*(1-y))
        f2=(xfa-x)/(y-x)
        f=[f1,f2]
        return f
    p=[x,y]
    p0=[0.3,0.7]
    p=fsolve(f,p0)    
    print p
    
    xr=p[0]
    yp=p[1]
    
    print xr
    print yp
    
    
    
    

if (k!=0 and k!=1):
    print ("wrong choice")

    
    
    







