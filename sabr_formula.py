# -*- coding: utf-8 -*-
"""
Created on Sat Oct  2 10:23:06 2021

@author: Apa
"""


import numpy as np
 

def SABR1_vols(K,Fwd, Texp,alpha,beta,nu,rho):
    ''' Returns Sabr volatilities for a given set of strikes. Based on original SABR formula.
        Wilmott magazine  - 2.17 and 2.18
    '''
    vols = np.zeros(len(K))
    
    for ii in range(0,len(K)):
        
        if K[ii] == Fwd: #ATM
            B1 = (1.0 - beta)**2.0*alpha**2.0/(24.0*Fwd**(2.0 - 2.0*beta)) #2.18
            B2 = rho*beta*alpha*nu/(4.0*Fwd**(1.0 - beta)) #2.18
            B3 = (2.0 - 3.0*rho**2)*nu**2.0/24.0 #2.18
            
            vols[ii] = (alpha/Fwd**(1 - beta))*(1 + (B1 + B2 + B3)*Texp )
        
        else:   #OTM
        
            logfK = np.log(Fwd/K[ii])
            fkbpow = (Fwd*K[ii])**((1.0 - beta)/2.0)
            z = nu*fkbpow*logfK/alpha  #2.17b
            xz = np.log((np.sqrt(1.0 - 2.0*rho*z + z**2.0 ) + z - rho)/(1.0-rho)) #2.17c
            
            A1 = ((1.0-beta)**2.0)*(alpha**2.0)/(24.0*fkbpow**2.0) #2.17a
            A2 = (rho*beta*nu*alpha)/(4.0*fkbpow) #2.17a
            A3 = (2.0-3.0*rho**2)*nu**2.0/24.0 #2.17a
            A4 = ((1.0-beta)**2)*(logfK**2) #2.17a
            A41 = A4/24.0
            A42 =  A4*A4/1920.0 #2.17a
            
            vols[ii] = (alpha*z*(1 + (A1 + A2 + A3)*Texp ))/(fkbpow*xz*(1 + A41 + A42 ))
            
    return vols