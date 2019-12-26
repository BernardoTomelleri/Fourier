# -*- coding: utf-8 -*-
"""
Created on Tue Dec 24 16:29:20 2019

@author: berna
"""
import numpy as np
from matplotlib import pyplot as plt
from collections import namedtuple
from scipy.optimize import curve_fit
from scipy import signal as sg

# Definizione componenti dell'integratore A
R_a = 6.72e3
dR_a = np.sqrt((R_a*8e-3)**2 + 1e-2**2) 
C_a = 2.2e-6
dC_a = C_a*0.1
tau_a = R_a*C_a
dtau_a = np.sqrt((R_a*dC_a)**2 + (C_a*dR_a)**2 +2*dR_a*dC_a)
fT_a = 1./(2*np.pi*tau_a)
wT_a = 2*np.pi*fT_a*1e-3
#dT_a = 1./2*np.pi*(dtau/tau_a**2)
dfT_a = fT_a*0.1
print('Frequenza Taglio A = %.f +- %.f' %(fT_a, dfT_a))
# Definizioni variabili per la serie di Fourier
t0 = 1         # tempi iniziale
tN = 99       # e finale
n_points = 1000 # punti equispaziati
sums = 50   # termini nella serie
# Array di appoggio per lo sviluppo in serie
tt=np.linspace(t0, tN, n_points)
c = np.zeros(sums)
w = np.zeros_like(c)
A = np.zeros_like(c)
Df = np.zeros_like(c)
# Tuple di appoggio per le funzioni
args=namedtuple('pars', 'f A_pp B phi')
funcs = []
# Parametri di default: frequenza = 1 [kHz], Ampiezza(p-p) = 1 [arb un],
# DC offset = 0 [arb un], fase = 0 [pi rad]
def sqw(t, f=1, A_pp=1, B=0, phi=0):
    fk = 0
    t+= phi/(2*f)
    for k in range(1, sums, 2):
        c[k] = 2./(k*np.pi)
        w[k] = k*2*np.pi*f
        fk += c[k]*np.sin(w[k]*t)
    t-= phi/(2*f)
    return fk*A_pp + B

def square(t, f=1, A_pp=1, B=0, phi=0, duty=0.5):
        return (A_pp/2.)*sg.square(np.pi*2*f*t + np.pi*phi, duty) + B

def trg(t, f=1, A_pp=1, B=0, phi=0):
    t+= phi/(2*f)
    fk = 0
    for k in range(1, sums, 2):
        c[k] = (2./(k*np.pi))**2
        w[k] = k*2*np.pi*f
        fk += c[k]*np.cos(w[k]*t)
    t-= phi/(2*f)
    return fk*A_pp + B

def triangle(t, f=1, A_pp=1, B=0, phi=0, duty=0.5):
        return (A_pp/2.)*sg.sawtooth(2*np.pi*f*t + np.pi+phi, duty) + B
    
def intg(f, fT=fT_a):
    return 1./np.sqrt(1.+(f/fT)**2) 
def fin(t, f=1, A_pp=1, B=0, phi=0):
    t+= phi/(2*f)
    fk = 0
    for k in range(1, sums, 2):
        c[k] = 2./(k*np.pi)
        w[k] = k*2*np.pi*f
        A[k] = intg(w[k], wT_a)
        Df[k] = np.arctan(-w[k]/wT_a)
        fk += A[k]*c[k]*np.sin(w[k]*t + Df[k])
    t-= phi/(2*f)
    return fk*A_pp + B    

def chitest(y, dy, model, ddof=0):
    res = y - model
    resnorm = res/dy
    ndof = len(y) - ddof
    chisq = (resnorm**2).sum()
    sigma = (chisq - ndof)/np.sqrt(2*ndof)
    return chisq, ndof, sigma    

# Variabili di controllo dei grafici in ordine decrescente di pesantezza
log=True
tex=True
tick=True
DSO = False
# Impostazioni opzionali per tipografia in LaTeX
if tex:
   plt.rc('text', usetex=True)
   plt.rc('font', family='serif')
dy = np.full(len(tt), 1)
init=(0.05, 734, 160, 0)
delta = sqw(tt, *init)-square(tt,*init)
chisq, ndof, sigma = chitest(sqw(tt, *init),dy, square(tt, *init))
print('Chi quadro/ndof = %f/%d [%+.1f]' % (chisq, ndof, sigma))
fig1,(ax1, ax2) = plt.subplots(2,1, True, gridspec_kw={'wspace':0.05,
     'hspace':0.05, 'height_ratios': [3, 1]})
ax1.grid(c = 'gray', ls = '--', alpha=0.7)
ax1.plot(tt, sqw(tt, *init), c='black',
        label='$\lambda = %d$ \n $\chi^2 = %.f/%d$' %(sums,chisq,ndof), zorder=10, alpha=0.7)
ax1.plot(tt, square(tt, *init), c='blue',
        label='scipy', zorder=10, alpha=0.7)
ax1.set_xlabel('Tempo $t$ [ms]', x=0.92)
ax1.set_ylabel('Ampiezza [arb. un]')
ax1.minorticks_on()
if tick:
    ax1.yaxis.set_major_locator(plt.MultipleLocator(100.))
    ax1.yaxis.set_minor_locator(plt.MultipleLocator(20))  
ax1.tick_params(direction='in', length=5, width=1., top=True, right=True)
ax1.tick_params(which='minor', direction='in', width=1., top=True, right=True)
legend = ax1.legend(loc ='best', framealpha = 0.3)
legend.set_zorder(100)

ax2.set_xlabel('Tempo $t$ [ms]', x=0.92)
ax2.set_ylabel('Residui')
ax2.axhline(0, c='r', lw=1.1, zorder=10)
ax2.errorbar(tt, delta, None, None, 'ko', elinewidth = 0.7,
              capsize=0.7, ms=1., ls='--', lw=0.8, zorder=0)
ax2.grid(color ='gray', ls = '--', alpha=0.7)
ax2.ticklabel_format(axis='both', style='sci', scilimits=None,
                     useMathText=True)
ax2.minorticks_on()
ax2.xaxis.set_major_locator(plt.MultipleLocator(10))
ax2.xaxis.set_minor_locator(plt.MultipleLocator(2))
if tick:
    ax2.yaxis.set_major_locator(plt.MultipleLocator(1e4/sums))
    ax2.yaxis.set_minor_locator(plt.MultipleLocator(2e3/sums))
ax2.tick_params(direction='in', length=5, width=1., top=True, right=True)
ax2.tick_params(which='minor', direction='in', width=1., top=True, right=True)

delta = trg(tt, *init)-triangle(tt,*init)
chisq, ndof, sigma = chitest(trg(tt, *init),dy, triangle(tt, *init))
print('Chi quadro/ndof = %f/%d [%+.1f]' % (chisq, ndof, sigma))
fig1,(ax1, ax2) = plt.subplots(2,1, True, gridspec_kw={'wspace':0.05,
     'hspace':0.05, 'height_ratios': [3, 1]})
ax1.grid(c = 'gray', ls = '--', alpha=0.7)
ax1.plot(tt, trg(tt, *init), c='black',
        label='$\lambda = %d$ \n $\chi^2 = %.f/%d$' %(sums,chisq,ndof), zorder=10, alpha=0.7)
ax1.plot(tt, triangle(tt, *init), c='blue',
        label='scipy', zorder=10, alpha=0.7)
ax1.set_xlabel('Tempo $t$ [ms]', x=0.92)
ax1.set_ylabel('Ampiezza [arb. un]')
ax1.minorticks_on()
if tick:
    ax1.yaxis.set_major_locator(plt.MultipleLocator(100.))
    ax1.yaxis.set_minor_locator(plt.MultipleLocator(20))  
ax1.tick_params(direction='in', length=5, width=1., top=True, right=True)
ax1.tick_params(which='minor', direction='in', width=1., top=True, right=True)
legend = ax1.legend(loc ='right', framealpha = 0.3)
legend.set_zorder(100)

ax2.set_xlabel('Tempo $t$ [ms]', x=0.92)
ax2.set_ylabel('Residui')
ax2.axhline(0, c='r', lw=1.1, zorder=10)
ax2.errorbar(tt, delta, None, None, 'ko', elinewidth = 0.7,
              capsize=0.7, ms=1., ls='--', lw=0.8, zorder=0)
ax2.grid(color ='gray', ls = '--', alpha=0.7)
ax2.ticklabel_format(axis='both', style='sci', scilimits=None,
                     useMathText=True)
ax2.minorticks_on()
ax2.xaxis.set_major_locator(plt.MultipleLocator(10))
ax2.xaxis.set_minor_locator(plt.MultipleLocator(2))
if tick:
    ax2.yaxis.set_major_locator(plt.MultipleLocator(1))
    ax2.yaxis.set_minor_locator(plt.MultipleLocator(0.2))
ax2.tick_params(direction='in', length=5, width=1., top=True, right=True)
ax2.tick_params(which='minor', direction='in', width=1., top=True, right=True)

t, V = np.loadtxt('../RC/data_int/syncsqw_50.txt', unpack=True)
F, dF, Att, dAtt, Btt, dBtt = np.genfromtxt('../RC/Att_data.csv', float, delimiter=',',
skip_header = 1, usecols=(0, 1, 8, 9, 10, 11), unpack = True)
if DSO:
    t, V = np.genfromtxt('./DSO1000/DSO00001.csv', float, delimiter=',',
                         skip_header = 2, usecols=(0, 2), unpack = True)

# Trasformazione dei dati nelle grandezze da fittare
x = t/1.e3
y = V
# Estrazione di un sottointervallo di dati
x_min = t0
x_max = tN
x1 = x[x>x_min]; sx = x1[x1<x_max];
y1 = y[x>x_min]; sy = y1[x1<x_max];
dx = np.full(len(sx), 4.)/1.e3
dy = np.full(len(sy), 2.)
# Fit dei campionamenti con una pinna di squalo
init=(0.05, 300, 160, 1)
popt, pcov = curve_fit(fin, sx, sy, init, dy, absolute_sigma = True)
f_fit, A_fit, B_fit, phi_fit = popt
df_fit, dA_fit, dB_fit, dphi_fit = np.sqrt(pcov.diagonal())
chisq, ndof, sigma = chitest(sy, dy, fin(sx, *popt), ddof=len(popt))
print('Chi quadro/ndof = %f/%d [%+.1f]' % (chisq, ndof, sigma))
fig,(ax1, ax2) = plt.subplots(2,1, True, gridspec_kw={'wspace':0.05,
     'hspace':0.05, 'height_ratios': [3, 1]})
ax1.errorbar(sx, sy, dy, dx,'ko', ms=1.5, elinewidth=1., capsize=1.5,
            ls='',label='data', zorder=5)
ax1.plot(tt, fin(tt, *popt),
        label='fin \n $\chi^2 = %.f/%d$' %(chisq, ndof), zorder=10, alpha=0.7)
ax1.set_xlabel('Tempo $t$ [ms]', x=0.92)
ax1.set_ylabel('Ampiezza [arb. un]')
ax1.grid(c = 'gray', ls = '--', alpha=0.7)
ax1.yaxis.set_major_locator(plt.MultipleLocator(20))
ax1.yaxis.set_minor_locator(plt.MultipleLocator(4))
ax1.tick_params(direction='in', length=5, width=1., top=True, right=True)
ax1.tick_params(which='minor', direction='in', width=1., top=True, right=True)
legend = ax1.legend(loc ='right', framealpha = 0.3)
res = sy - fin(sx, *popt)
resnorm = res/dy
ax2.set_xlabel('Tempo $t$ [ms]', x=0.92)
ax2.set_ylabel('Residui')
ax2.axhline(0, c='r', zorder=10)
ax2.errorbar(sx, resnorm, None, None, 'ko', elinewidth = 0.7,
              capsize=0.7, ms=1., ls='--', lw=0.8, zorder=0)
ax2.grid(color ='gray', ls = '--', alpha=0.7)
ax2.ticklabel_format(axis='both', style='sci', scilimits=None,
                     useMathText=True)
ax2.xaxis.set_major_locator(plt.MultipleLocator(5))
ax2.xaxis.set_minor_locator(plt.MultipleLocator(1))
ax2.yaxis.set_major_locator(plt.MultipleLocator(2))
ax2.yaxis.set_minor_locator(plt.MultipleLocator(0.4))
ax2.tick_params(direction='in', length=5, width=1., top=True, right=True)
ax2.tick_params(which='minor', direction='in', width=1., top=True, right=True)
print('f = %f +- %f Hz' % (f_fit*1e3, df_fit*1e3))
print('A = %f +- %f mV' % (A_fit, dA_fit))
print('B = %f +- %f mV' % (B_fit, dB_fit))
print('phi = %f +- %f rad' % (phi_fit, dphi_fit))
plt.show()