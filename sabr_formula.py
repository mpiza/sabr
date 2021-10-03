# -*- coding: utf-8 -*-
"""
Created on Sat Oct  2 10:23:06 2021

@author: Apa
"""


import numpy as np
import math

def SABR_normal_vols(K,Fwd, Texp,alpha,beta,nu,rho):
    
    
   beta_bar = 1 - beta
   mny = Fwd - K
   mny2 = np.multiply(mny, mny)
   fk_05b = np.power(Fwd * K, -0.5 * beta)
   z = nu * mny / alpha * fk_05b
   z2 = np.multiply(z,z)
   sqrt_z2_zr_1 = np.sqrt(z2 - 2 * rho *z + 1)
   xz = np.log((sqrt_z2_zr_1 + z - rho) / (1 - rho))
   fb_kb = math.pow(Fwd, beta_bar) - np.power(K, beta_bar)
   temp = np.divide(beta_bar * nu * mny2, fb_kb)
   temp = np.multiply(temp, fk_05b)
   scalar = np.divide(temp, xz)

   fk_05b_bar = np.power(Fwd*K, -0.5 * beta_bar)
   fk_b_bar = np.multiply(fk_05b_bar,fk_05b_bar)
   
   term1 = beta * (beta - 2) * alpha * alpha / 24.0 * fk_b_bar
   term2 = alpha * beta * rho * nu * 0.25 * fk_05b_bar
   term3 = (2 - 3 * rho * rho) / 24.0 * nu * nu
    
   for ii in range(0,len(K)):
        
        if K[ii] == Fwd: #ATM
        
        
        
            beta_bar = 1 - beta
            k_beta_bar = np.power(K, -beta_bar)
            term1 = beta * (beta - 2) * alpha * \
            alpha / 24.0 * np.multiply(k_beta_bar, k_beta_bar)
            term2 = alpha * beta * rho * nu * 0.25 * k_beta_bar
            term3 = (2 - 3 * rho * rho) / 24.0 * nu * nu
            scalars = alpha * np.power(K, beta)
        
            vols = np.multiply(scalars, 1 + (term1 + term2 + term3) * Texp)
            
             
        
        else:   #OTM
        
           vols = np.multiply(scalar, 1 + (term1 + term2 + term3) * Texp)
            
        return vols   
           
def SABR_normal_vols_shifted(K,Fwd, Texp,alpha,beta,nu,rho,shift):
    
     vols =  SABR_normal_vols(K + shift,Fwd + shift, Texp,alpha,beta,nu,rho)
     
     print(vols)
     
     return vols

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

def SABR1_vols_shifted(K,Fwd, Texp,alpha,beta,nu,rho,shift):
    ''' Returns Sabr volatilities for a given set of strikes. Based on original SABR formula.
        Wilmott magazine  - 2.17 and 2.18
    '''
    vols = SABR1_vols(K + shift,Fwd + shift, Texp,alpha,beta,nu,rho) 
    return vols