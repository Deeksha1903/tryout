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
Pp=19.0 # permeate pressure in cmHg
Pf=190.0 # feed pressure in cmHg
r=Pf/Pp # pressure ratio
xfa=0.21 # mole fraction of oxygen in the feed
L=2.54*10**-3 # lenght of the membrane

np=theta*nf
nr=nf-np

# we assume completely mixed model to get the initial guess values to solvethe differental equations


def equation(yp):
    f=yp/(yp-1)-(a*(xfa-r*yp))/((1-xfa)-r*(1-yp))
    return f
A=fsolve(equation,0.3)

yp=A[0]  # initial permeate composition at x=xfa
print "initial permeate composition"
print yp

npa=theta*nf*yp # rate of transger of component A
Na=(Pma/L)*(xfa*Pf-yp*Pp) # flux of transfer of A
print "Flux through the membrane is"
print Na
Am=npa/Na # area of the membrane
Am1=Am*10**-4 # area required in m^2
print "Area of membrane in m^2 is "
print Am1
Pmb1=Pmb/L
Kf=(nf/Pmb1*Am*Pp)*(10**-6*10**-4*10**-4/10**-6*10**-4*10**-2)


k=input("0-cocurrent 1-countercurrent")

if (k==0):
    " for cocurrent flow "
    
    
    def deriv(z,A):
        [x,y]=z
        
        dxdA=(1/Kf)*((x-y)/(y-xfa))*(a*(1-x)*(x*r-y)-x*((1-x)*r)-(1-y))
        dydA=(1/Kf)*((x-y)/(x-xfa))*(a*(1-y)*(x*r-y)-y*((1-x)*r)-(1-y))
        
        ''' reference: A simple analysis for gas separation membranes by 
        Richard A. Davis and Orville E. Sandall
        university of Minnesota Duluth'''
        return [dxdA,dydA]
    
    # initial conditions
    x0=xfa+0.00001
    y0=yp
    z0=[x0,y0]
    print z0
    # time grid for integration
    A=scipy.linspace(0.0,1.0,11.0)
    
    
    p=scipy.integrate.odeint(deriv, z0, A)
    
    x=p[:,0]
    y=p[:,1]
    soln=scipy.array([[A],[x],[y]])
    
    print "COCURRENT FLOW MODEL"
    print soln
    
    plt.plot(A,x,'g')
    plt.xlabel('length')
    plt.ylabel('retentate composition')
   
    plt.show()

    plt.plot(A,y,'b')
    plt.xlabel('length')
    plt.ylabel('permeate composition')
    
    plt.show()
    
    print "COMMENTS: We start iterating from left to right of the membrane."
    print "The feed and permeate enter at the left end of the membrane"
    print "Across the length the permeate composition increases and retentate composition decreases ; shown by the graph"
    print "Desired sepaartion is achieved in cocurrent flow model"
    

if (k==1):
    "for countercurrent flow "
    
    
    xr=0.208 #desired retentate composition
    print "retentate composition which is desired"
    print xr
    
    def equation(yg):
        f=yg/(yg-1)-(a*(xr-r*yg))/((1-xr)-r*(1-yg))
        return f
    A=fsolve(equation,0.3)
    
  
    yg=A[0]
    print "permeate composition at the right entrance"
    print yg  # initial permeate composition at x=xr
    
    
    Kr=(nr/Pmb1*Am*Pp)*(10**-6*10**-4*10**-4/10**-6*10**-4*10**-2)
    
    
    def deriv(z,A):
        [x,y,n]=z
        
        dxdA=(1/Kr)*((x-y)/(xr-yg))*(a*(1-x)*(x*r-y)-x*((1-x)*r)-(1-y))
        dydA=(1/Kr)*((x-y)/(xr-x))*(a*(1-y)*(x*r-y)-y*((1-x)*r)-(1-y))
        dndA=(1/Kr)*(a*(x*r-y)+(1-x)*r-(1-y))
        
        ''' reference: A simple analysis for gas separation membranes by 
        Richard A. Davis and Orville E. Sandall
        university of Minnesota Duluth'''
        return [dxdA,dydA,dndA]
    
    # initial conditions
    x0=xr-0.001
    y0=yg
    n0=1.0
    
    z0=[x0,y0,n0]
    print z0
    # time grid for integration
    A=scipy.linspace(0.0,1.0,11.0)
    
    
    p=scipy.integrate.odeint(deriv, z0, A)
    
    x=p[:,0]
    y=p[:,1]
    n=p[:,2]
    
    soln=scipy.array([[A],[x],[y]])
    
    print "COUNTERCURRENT FLOW MODEL"
    print soln
    
    plt.plot(A,x,'g')
    plt.xlabel('length')
    plt.ylabel('retentate composition')
    
    plt.show()

    plt.plot(A,y,'b')
    plt.xlabel('length')
    plt.ylabel('permeate composition')
    
    plt.show()
    
    print "COMMENTS: We start iterating from right to left of the membrane."
    print "The feed enters at the left and permeate enters at the right end of the membrane"
    print "Across the length the permeate composition  and retentate composition increase ; shown by the graph"
    print "Desired sepaartion is not achieved in countercurrent flow model"
    
    
        
    
if (k!=0 and k!=1):
    print ("wrong choice")

    
    
    







