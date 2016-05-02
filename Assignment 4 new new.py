# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 21:04:58 2016

@author: Kaivalya
"""

import scipy
import matplotlib.pyplot 


class MassTransfer:
    
    #xl= win32com.client.gencache.EnsureDispatch("Excel.Application")
    #wb=xl.Workbooks('Data.xlsx').............work book Data,xlsx
    #sheet=wb.Sheets('Data')..................sheet Data
    
    #def getdata(sheet, Range):
        #data= sheet.Range(Range).Value
        #data=scipy.array(data)
        #data=data.reshape((1,len(data)))[0]
        #return data
#=================================================================================================#    
    '''Input Values'''    
    Lwin=1.0 #kmoles  liquid flow rate of water IN
    Gmin=0.5 #kmoles  gas glow rate of methane IN
    Gcin=0.5 #kmoles  gas flow rate of CO2 IN
    Lcin=0.0 #kmoles  Liquid flow rate of CO2 IN
    Lmin=0.0 #kmoles  liquid flow rate of methane IN
    Klac=18 # m/hr..............................Ref 5 
    Klam=102.9 # m/hr...........................Ref 5
    Kgaw=0.001367*3600 #m/hr
    Hc=3.35*(10**-4) #mol/m^3.Pa..................Ref 4
    Hm=1.38*(10**-4) #mol/m^3.Pa..................Ref 4
        
    def psat(T):
        Psat=((10**(7.9668-(1668.21/(T+228.01))))/760)*100000
        return Psat
    
    Pop=1.013*(10**5)
    Top=298.16
    R=8.314
    A=0.05
    L=1
#=================================================================================================#
    '''First set of guess values'''
    Gcout1=0.30
    Gmout1=0.30
    Gwout1=0.40
#=================================================================================================#
    '''Second set of guess values'''
    Gcout2=0.51
    Gmout2=0.51
    Gwout2=0.01
#=================================================================================================#
    '''Initial values of moles in liquid stream for each component for both set of guess values'''
    Lm1=Lmin
    Lc1=Lcin
    Lw1=Lwin
    
    Lm2=Lmin
    Lc2=Lcin
    Lw2=Lwin
#=================================================================================================#
    Gmout2new=Gmout2
    Gcout2new=Gcout2
    Gwout2new=Gwout2
    n=10
#=================================================================================================#
    '''Creating a list of different N values'''
    list_n=scipy.array([10*i for i in range(1,10)])
    calcerr=[]
    for n in list_n:
#=================================================================================================#
        dL=L/(n-1)
        Err2=1
        Err=1
#=================================================================================================#  
        '''Initiating the while loop'''
        while abs(Err2)>0.001:
            Gm1=Gmout1
            Gc1=Gcout1
            Gw1=Gwout1
            
            Gm2=Gmout2new
            Gc2=Gcout2new
            Gw2=Gwout2new
            for i in range(n):
                delGm1=Klam*A*(-Pop*(Gm1/(Gm1+Gc1+Gw1))*Hm+(Lm1/(Lw1*18*0.001)))
                delGc1=Klac*A*(-Pop*(Gc1/(Gm1+Gc1+Gw1))*Hc+(Lc1/(Lw1*18*0.001)))
                delGw1=Kgaw*A*(psat(25)-Pop*(Gw1/(Gm1+Gc1+Gw1)))
                Gm1new=Gm1-delGm1*dL
                Lm1new=Lm1+delGm1*dL
                Gc1new=Gc1-delGc1*dL
                Lc1new=Lc1+delGc1*dL
                Gw1new=Gw1+delGw1*dL
                Lw1new=Lw1-delGw1*dL
                
                delGm2=Klam*A*(Pop*(Gm2/(Gm2+Gc2+Gw2))*Hm-(Lm2/(Lw2*18*0.001)))
                delGc2=Klac*A*(Pop*(Gc2/(Gm2+Gc2+Gw2))*Hc-(Lc2/(Lw2*18*0.001)))
                delGw2=Kgaw*A*(psat(25)-Pop*(Gw2/(Gm2+Gc2+Gw2)))
                Gm2new=Gm2-delGm2*dL
                Lm2new=Lm2+delGm2*dL
                Gc2new=Gc2-delGc2*dL
                Lc2new=Lc2+delGc2*dL
                Gw2new=Gw2+delGw2*dL
                Lw2new=Lw2-delGw2*dL
                
                Gm1=Gm1new
                Gc1=Gc1new
                Gw1=Gw1new
                Lm1=Lm1new
                Lc1=Lc1new
                Lw1=Lw1new
                
                Gm2=Gm2new
                Gc2=Gc2new
                Gw2=Gw2new
                Lm2=Lm2new
                Lc2=Lc2new
                Lw2=Lw2new
#=================================================================================================#
            '''Error calculation for first set of guess values'''    
            Em1=(0.5-Gm1)
            Ec1=(0.5-Gc1)
            Ew1=(0.0-Gw1)
            E1=[Em1,Ec1,Ew1]
            #print(E1)
            Err1=Em1+Ec1+Ew1
            #print(Err1)
            
            '''Error calculation for second set of guess values'''
            Em2=(0.5-Gm2)
            Ec2=(0.5-Gc2)
            Ew2=(0.0-Gw2)
            E2=[Em2,Ec2,Ew2]
            #print(E2)
            Err2=Em2+Ec2+Ew2
            #print(Err2)
#================================================================================================#
            '''Modification of second set of guess values by Secant Method'''        
            Gmout2new=Gmout2+((Gmout2-Gmout1)/(Em2-Em1))*Em2
            Gcout2new=Gcout2+((Gcout2-Gcout1)/(Ec2-Ec1))*Ec2
            Gwout2new=Gwout2+((Gwout2-Gwout1)/(Ew2-Ew1))*Ew2
        Massin=16*(0+0.5)+44*(0+0.5)+18*(1+0)
        Massout=16*(Gmout2+Lm2)+44*(Gcout2+Lc2)+18*(Gwout2+Lw2)
        Err=abs(Massin-Massout)/Massin
        calcerr+=[abs(Err2)]
#================================================================================================#
    print (list_n)
    print (calcerr)
    '''For plotting graph'''
    matplotlib.pyplot.plot(list_n,calcerr)
    matplotlib.pyplot.show()
#=================================================================================================#

'''
References:
1.Perrys handbook 1999 page no 2114,2115
2.Korean J. Chem. Eng., 32(6), 1060-1063 (2015) 
3.Convective Mass trasnfer Chapter 3
4.Henry-3.0 pdf,R.Sander: Henrys law constant
5.International Journal of Scientific and Research Publications, Volume 4, Issue 4, April 2014.
'''