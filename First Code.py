from pylab import *
import scipy
from scipy.integrate import odeint
from scipy.optimize import fsolve
import pylab
import numpy

kla1=(10**-2.6)/17.226
kla2=(0.2/3600)
a=5.0
h1=(1.0/0.00034)
h2=(1.0/0.000014)
psat=0.0231*101325.0
z=5.0
n=100.0
def contact(f,z):
    
    c1=kla1*a*(p/h1*(f[0]/(f[0]+f[1]+f[2])-f[3]*1000/18/f[5]))
    c2=kla2*a*(p/h2*(f[1]/(f[0]+f[1]+f[2])-f[4]*1000/18/f[5]))
    c3=kga*a*(p*f[2]/(f[0]+f[1]+f[2])-psat)
    return [-c1,-c2,-c3,-c1,-c2,-c3]
    
def shooting(l):
    
    z=numpy.linspace(0,10,101)
    f=odeint(contact,([50,50,0,l[0],l[1],l[2]]),z)
    e=[(f[100])[3],(f[100])[4],(f[100])[5]-100]
    return e
l=fsolve(shooting,[0.1169,0.001966,99.73013])
z=numpy.linspace(0,10,101)
f=odeint(contact,([50,50,0,l[0],l[1],l[2]]),z)
b=[[(f[100])[0],(f[100])[1],(f[100])[2],(f[100])[3],(f[100])[4],(f[100])[5]],
    [(f[0])[0],(f[0])[1],(f[0])[2],(f[0])[3],(f[0])[4],(f[0])[5]]]
print b
#pylab.plot(z,f[:,0])
#pylab.show()
#pylab.plot(z,f[:,1])
#pylab.show()
#pylab.plot(z,f[:,2])
#pylab.show()
#pylab.plot(z,f[:,3])
#pylab.show()
#pylab.plot(z,f[:,4])
#pylab.show()
#pylab.plot(z,f[:,5])
#pylab.show()
