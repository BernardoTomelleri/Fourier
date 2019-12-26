# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 12:57:47 2019

@author: berna
"""
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from collections import namedtuple
import sys
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
# Definizione componenti del derivatore B
R_b = 67.5
dR_b = np.sqrt((R_b*8e-3)**2 + 0.1**2) 
C_b = 0.2e-6
dC_b = C_b*0.1
tau_b = R_b*C_b
dtau_b = np.sqrt((R_b*dC_b)**2 + (C_b*dR_b)**2 +2*dR_b*dC_b)
fT_b = 1./(2*np.pi*tau_b)
wT_b = 2*np.pi*fT_b*1.e-3
#dT_b = 1./2*np.pi*(dtau/tau_a**2)
dfT_b = fT_b*0.1
print('Frequenza Taglio B = %.e +- %.e' %(fT_b, dfT_b))
# Definizioni variabili per la serie di Fourier
t0 = 0.         # tempi iniziale
tN = 100.       # e finale
n_points = 1000 # punti equispaziati
sums = 10000    # termini nella serie
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

def trg(t, f=1, A_pp=1, B=0, phi=0):
    t+= phi/(2*f)
    fk = 0
    for k in range(1, sums, 2):
        c[k] = (2./(k*np.pi))**2
        w[k] = k*2*np.pi*f
        fk += c[k]*np.cos(w[k]*t)
    t-= phi/(2*f)
    return fk*A_pp + B

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

def derv(f, fT_int=fT_a, fT_der=fT_b):
    return intg(f, fT_int)*1./np.sqrt(1.+(fT_der/f)**2)
def cas(t, f=1, A_pp=1, B=0, phi=0):
    t+= phi/(2*f)
    fk = 0
    for k in range(1, sums, 2):
        c[k] = 2./(k*np.pi)
        w[k] = k*2*np.pi*f
        A[k] = derv(w[k], wT_a, wT_b)
        Df[k] = np.arctan(-w[k]/wT_a) + np.arctan(wT_b/w[k])
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
log=False
tex=True
tick=False
DSO = False
# Impostazioni opzionali per tipografia in LaTeX
if tex:
   plt.rc('text', usetex=True)
   plt.rc('font', family='serif')
   
def plot(waves, x, args=args(f=1.,A_pp=1., B=0., phi=0.), tick=False):
    # Grafico di segnali nel tempo
    fig, ax = plt.gcf(), plt.gca()
    for f in waves:
        ax.plot(x, f(x, *args))
    ax.grid(c = 'gray', ls = '--', alpha=0.7)
    ax.set_xlabel('Tempo $t$ [ms]', x=0.92)
    ax.set_ylabel('Ampiezza [arb. un]')
    ax.minorticks_on()
    if tick:
        ax.xaxis.set_major_locator(plt.MultipleLocator((tN-t0)/10.))
        ax.xaxis.set_minor_locator(plt.MultipleLocator((tN-t0)/50.))
    ax.tick_params(direction='in', length=5, width=1., top=True, right=True)
    ax.tick_params(which='minor', direction='in', width=1., top=True, right=True)

# Extrazione dei vettori di campionamenti ed attenuazioni misurate
t, V = np.loadtxt('../RC/data_int/syncsqw_70.txt', unpack=True)
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
init=(0.07, 734, 160, 1)
popt, pcov = curve_fit(fin, sx, sy, init, dy, absolute_sigma = True)
f_fit, A_fit, B_fit, phi_fit = popt
df_fit, dA_fit, dB_fit, dphi_fit = np.sqrt(pcov.diagonal())
chisq, ndof, sigma = chitest(sy, dy, fin(sx, *popt), ddof=len(popt))
print('Chi quadro/ndof = %f/%d [%+.1f]' % (chisq, ndof, sigma))
fig, ax = plt.subplots()
ax.errorbar(sx, sy, dy, dx,'ko', ms=1.5, elinewidth=1., capsize=1.5,
            ls='',label='data', zorder=5)    
plot([fin], tt, args(*popt), tick=True)
# Attenuazione dell'onda simulata su un intervallo di frequenze frng
def sim(func=fin, x=sx, frng=F):
    simA=np.zeros(len(frng))
    for i in range(len(frng)):
        simA[i] = (max(func(x, frng[i]*1.e-3))- min(func(x, frng[i]*1.e-3)))
    return simA
if not log:
    plt.show()
    sys.exit()
fig, ax = plt.subplots()
ax.set_yscale('log')
ax.set_xscale('log')
FF = np.logspace(np.log10(min(F)), np.log10(max(F))+0.2, 400)
A_sim, AB_sim = sim(fin, x=tt, frng=FF), sim(cas, tt, FF)
xx = np.logspace(np.log10(min(F)), np.log10(max(F))+0.2, 2000)
ax.errorbar(F, Att, dAtt, dF,'ko', ms=1.5, elinewidth=1., capsize=1.5,
            ls='',label='data$_A$', zorder=10)
ax.errorbar(F, Btt, dBtt, dF,'bx', ms=1.5, elinewidth=1., capsize=1.5,
            ls='',label='data$_{A+B}$', zorder=10)
cA, nA, sA = chitest(Att, dAtt, model=sim(fin, x=sx, frng=F))
ax.errorbar(FF, A_sim, None, None,'r^', ms=0.5, ls='-', lw=1.,
            label='sim$_A$\n $\chi^2 = %.f/%d$' %(cA,nA), zorder=2)
cB, nB, sB = chitest(Btt, dBtt, model=sim(cas, x=sx, frng=F))
ax.errorbar(FF, AB_sim, None, None,'gv', ms=0.5, ls='-', lw=1.,
            label='sim$_{A+B}$\n $\chi^2 = %.f/%d$' %(cB,nB), zorder=2)
ax.axhline(fT_a/fT_b, c='gray', ls='-.',lw=1.2)
ax.annotate('$f_{TA}/f_{TB}$',(6e4, fT_a/fT_b +4e-4))
ax.plot(xx, intg(xx, fT_a), c='gray', ls='--', lw=1.2, zorder=1)
ax.annotate('$1/\sqrt{1+(f/f_{T_A})^2}$',(3.5, intg(10)-0.61))
ax.plot(xx, derv(xx, fT_a, fT_b), c='gray', ls='--', lw=1.2, zorder=1)    
ax.annotate('$1/{\sqrt{1+(f/f_{T_A})^2} \sqrt{1+(f_{T_B}/f)^2}}$',
            (7, derv(100)-7e-4))
ax.set_ylabel('Gain $A(f)$', y=0.9)
ax.set_xlabel('Frequency $f$ [Hz]', x=0.9)
ax.set_title('$A$: LPF with $f_{T_A} = 10$ [Hz], $B$: HPF with $f_{T_B} = 10$ [kHz]')
ax.tick_params(direction='in', length=5, width=1., top=True, right=True)
ax.tick_params(which='minor', direction='in', width=1., top=True, right=True)
legend = ax.legend(loc ='best', framealpha = 0.4)
plt.show()