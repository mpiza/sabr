# -*- coding: utf-8 -*-
"""
Created on Sat Oct  2 10:49:18 2021

@author: Apa
"""

import sabr_formula 
import numpy as np
from scipy.optimize import curve_fit

def atm_sigma_to_alpha(Fwd,Texp,sigma_atm,beta,nu,rho):
    '''Returns alpha given forward, ATM volatility,expiry (Texp),
    beta, nu and rho.  Solves  cubic equation for alpha, equation (2.18) 
    "Managing smile risk." The Best of Wilmott 1 (2002): 249-296. 
    '''
    #Coeffceints of 3rd order polynomialin alpha (2.18)
    p_3 = -sigma_atm
    p_2 =  (1 + (2-3*rho**2)*nu**2*Texp/24)/Fwd**(1.-beta)
    p_1 = rho*beta*nu*Texp/(4*Fwd**(2-2*beta))
    p_0 = (1-beta)**2*Texp/(24*Fwd**(3-3*beta))
    coeffs = [p_0,p_1,p_2,p_3]
    
    r = np.roots(coeffs)    #find the roots of the cubic equation
    
    return r[(r.imag==0) & (r.real>=0)].real.min() 

def SABR_calibration_fix_atm(Fwd, Texp, sigma_atm, beta, strikes, vols,guess):
    ''' Returns the parameters alpha, nu and rho given beta, 
    forward, a list of market volatilities and  strike 
    spreads. Instead of doing optimization in all three parameters, 
    calculate alpha from nu and rho and optimize in rho and nu.
    '''
    def func_to_optimize(K,nu,rho):
        alpha = atm_sigma_to_alpha(Fwd,Texp,sigma_atm,beta,nu,rho)
        return  sabr_formula.SABR1_vols(K,Fwd,Texp,alpha,beta,nu,rho)
     
    method = 'lm'  #'trf'    
    opt_parameters, pcov = curve_fit(func_to_optimize, strikes, vols, p0 = (guess[1],guess[2]), maxfev=10000, method = method)
    print(pcov)  
    
    
    nu = opt_parameters[0] 
    rho = opt_parameters[1]
    alpha = atm_sigma_to_alpha(Fwd,Texp,sigma_atm,beta,nu,rho)
 
    res = sabr_formula.SABR1_vols(strikes,Fwd,Texp,alpha,beta,nu,rho) - vols
    print(res)
    res2 = np.sum(res**2)/res.size
    
    
#    nu = 0.00001;
#    rho = 0.000001;
#    alpha = 0.3
    
    return [alpha, nu, rho, res2]