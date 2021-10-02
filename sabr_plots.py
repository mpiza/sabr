# -*- coding: utf-8 -*-
"""
Created on Sat Oct  2 10:58:48 2021

@author: Apa
"""

import sabr_formula
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

Fwd = 0.03
Texp = 1.0

strikes_in_bps = np.linspace(-200,200,60)
strikes = Fwd + strikes_in_bps*0.0001


alpha0 = 0.2
beta0 = 1
nu0 = 0.1
rho0 = -0.5


vols = sabr_formula.SABR1_vols(strikes,Fwd,Texp,alpha0,beta0,nu0,rho0)

 

fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.25)
ax.set_ylim(0, 3*alpha0*Fwd**(1-beta0))


l, = plt.plot(strikes_in_bps, vols)
ax.margins(x=0)

axcolor = 'lightgoldenrodyellow'
ax_alpha = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
ax_beta = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
ax_nu = plt.axes([0.25, 0.01, 0.65, 0.03], facecolor=axcolor)
ax_rho = plt.axes([0.25, 0.05, 0.65, 0.03], facecolor=axcolor)

salpha = Slider(ax_alpha, 'alpha', 0.1, 1, valinit= alpha0)
sbeta = Slider(ax_beta, 'beta', 0.1, 1.0, valinit= beta0)
snu = Slider(ax_nu, 'nu', 0.001, 3, valinit= nu0)
srho = Slider(ax_rho, 'rho', -0.99, 1.0, valinit= rho0)

def update(val):
    alpha = salpha.val
    beta = sbeta.val
    nu = snu.val
    rho = srho.val
    ax.set_ylim(0, 3*alpha*Fwd**(1-beta))
    
    l.set_ydata(sabr_formula.SABR1_vols(strikes,Fwd,Texp,alpha,beta,nu,rho))
    fig.canvas.draw_idle()

salpha.on_changed(update)
sbeta.on_changed(update)
snu.on_changed(update)
srho.on_changed(update)

resetax = plt.axes([0.01, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')


def reset(event):
    salpha.reset()
    sbeta.reset()
    snu.reset()
    srho.reset()
    
    
button.on_clicked(reset)

#rax = plt.axes([0.025, 0.5, 0.15, 0.15], facecolor=axcolor)
#radio = RadioButtons(rax, ('red', 'blue', 'green'), active=0)


def colorfunc(label):
    l.set_color(label)
    fig.canvas.draw_idle()
#radio.on_clicked(colorfunc)


plt.show()