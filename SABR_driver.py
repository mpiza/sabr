# -*- coding: utf-8 -*-
"""
@author: William John Trenberth
email: w.j.trenberth@gmail.com
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sabr_formula
import sabr_calibration


def main(): 
     
    beta = 0.5
    Fwd =  0.04
    Texp = 5
    tenor = 1
    
    #A list of market volatilities at strikes corropsponding to strikes_in_bps below. 
    #sigmas = np.array([0.4040, 0.3541, 0.3218, 0.3107, 0.3048, 0.2975, 0.2923, 0.2873, 0.2870])
    sigmas = np.array([0.4,0.35, 0.33, 0.32, 0.3, 0.1, 0.3, 0.4, 0.5, 0.6,0.7])
    #The 'At the money volatility', corrosponding to a strike equal to the current forward price.
    #atm_sigma = 0.3048
    atm_sigma = 0.1
    #A list of strikes in bps (=0.0001) corrosponding to volatilites in sigmas
    strikes_in_bps = np.array([-200, -150,-100,-50,-25,0,25,50,100,150,200])
    #An inital guess of the parameters alpha, nu and rho.
    guess = [0.4, 3, 0.6]
    strikes = Fwd + strikes_in_bps*0.0001
    
     
    alpha, nu, rho, res2 = sabr_calibration.SABR_calibration_fix_atm(Fwd, Texp, atm_sigma, beta, strikes, sigmas, guess)
    print ("mse:", np.sqrt(res2))
    
    #PLOT
    Ks_in_bps = np.linspace(-200,200,60)
    Ks = Fwd + Ks_in_bps*0.0001
    vols_from_Ks = sabr_formula.SABR1_vols(Ks,Fwd,Texp,alpha,beta,nu,rho)
    textbox = "\n".join((r"$\alpha=$"+f"{round(alpha,6)}",r"$\beta=$"+f"{beta}",
                        r"$\rho=$"+f"{round(rho,6)}", r"$\nu=$"+f"{round(nu,6)}"))
    fig, ax = plt.subplots()
    plt.plot(strikes_in_bps, sigmas, 'x')
    plt.plot(Ks_in_bps,vols_from_Ks)
    plt.xlabel("Strikes in bps")
    plt.ylabel("Market volatilities")
    plt.title(f"{Texp} year into {tenor} year swaption")
    plt.text(0.6, 0.9, textbox, transform=ax.transAxes, fontsize=10,
        verticalalignment='top',bbox=dict(facecolor='white', alpha=0.7))
    
    #Saving the plot if desired.
    #plt.savefig(f"{t_exp} year into {tenor} year swaption"+".png")        



if __name__ == "__main__": main()

    
    
    
    