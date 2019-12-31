# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 19:52:17 2019

@author: berna
"""
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
# Definizioni variabili per la serie di Fourier
t0 = 0.0002        # tempi iniziale
tN = 0.08       # e finale
n_points = 1000 # punti equispaziati
sums = 1000    # termini nella serie

wT = 0.067 # frequenza angolare di taglio

# Array di appoggio per lo sviluppo in serie
tt=np.linspace(t0, tN, n_points)
c = np.zeros(len(tt))
w = np.zeros_like(c)
A = np.zeros_like(c)
Df = np.zeros_like(c)

def cas(t, f=1, A_pp=1, B=0, phi=0, wT_a=1, wT_b=1):
    t+= phi/(2*f)
    fk = 0
    for k in range(1, sums, 2):
        c[k] = 2./(k*np.pi)
        w[k] = k*2*np.pi*f
        A[k] = 1./(np.sqrt(1 + (w[k]/wT_a)**2)*np.sqrt(1 +(wT_b/w[k])**2))
        Df[k] = np.arctan(-w[k]/wT_a) + np.arctan(wT_b/w[k])
        fk += A[k]*c[k]*np.sin(w[k]*t + Df[k])
    t-= phi/(2*f)
    return fk*A_pp + B

def cas(t, f=1, A_pp=1, B=0, phi=0, wT_a=0.067, wT_b=74):
    t+= phi/(2*f)
    fk = 0
    for k in range(1, sums, 2):
        c[k] = 2./(k*np.pi)
        w[k] = k*2*np.pi*f
        A[k] = 1./(np.sqrt(1 + (w[k]/wT_a)**2)*np.sqrt(1 +(wT_b/w[k])**2))
        Df[k] = np.arctan(-w[k]/wT_a) + np.arctan(wT_b/w[k])
        fk += A[k]*c[k]*np.sin(w[k]*t + Df[k])
    t-= phi/(2*f)
    return fk*A_pp + B

t, DV = np.genfromtxt('./DSO/DSO00001.csv', float, delimiter=',',
                     skip_header = 2, usecols=(0, 1), unpack = True)

""" Extrazione dei dati corrispondenti all'intervallo di ddp selezionato
tempo in millisecondi, differenza di potenziale in unità arbitrarie """
t_min = t0
t_max = tN
t1 = t[t>t_min]; st = t1[t1<t_max];
DV1 = DV[t>t_min]; sDV = DV1[t1<t_max];
dt = np.full(len(st), 1e-5)
dDV = np.full(len(sDV), 9e-2)

log=False
tick=True
tex=True

init=(40.6, 3735, 0., 2., 19421, 54.7)

pars, covm = curve_fit(cas, st, sDV, init, dDV, absolute_sigma = False)
f_fit, A_fit, B_fit, phi_fit, wTa_fit, wTb_fit = pars
df_fit, dA_fit, dB_fit, dphi_fit, dwTa_fit, dwTb_fit = np.sqrt(covm.diagonal())
print('Parametri del fit:\n', pars)
print('Matrice di Covarianza:\n', covm)
print('f = %f +- %f Hz' % (f_fit, df_fit))
print('A = %f +- %f mV' % (A_fit, dA_fit))
print('B = %f +- %f mV' % (B_fit, dB_fit))
print('phi = %f +- %f rad' % (phi_fit, dphi_fit))
print('wTa = %f +- %f rad/s' % (wTa_fit, dwTa_fit))
print('wTb = %f +- %f rad/s' % (wTb_fit, dwTb_fit))
#Test Chi quadro per cas
res = sDV - cas(st, *pars)
resnorm = res/dDV
ndof = len(sDV) - len(init)
chisq = (resnorm**2).sum()
chirid = chisq/ndof
sigma = (chisq - ndof)/np.sqrt(2*ndof)
print('Chi quadro/ndof = %f/%d [%+.1f]' % (chisq, ndof, sigma))
print('Chi quadro ridotto_v:', chirid)

# Covarianza tra i parametri di cas
corr_fA = covm[0][1]/(df_fit*dA_fit)
corr_fB = covm[0][2]/(df_fit*dB_fit)
corr_fphi = covm[0][3]/(df_fit*dphi_fit)
corr_fwTa = covm[0][4]/(df_fit*dwTa_fit)
corr_fwTb = covm[0][5]/(df_fit*dwTb_fit)
corr_AB = covm[1][2]/(dA_fit*dB_fit)
corr_Aphi = covm[1][3]/(dA_fit*dphi_fit)
corr_AwTa = covm[1][4]/(dA_fit*dwTa_fit)
corr_AwTb = covm[1][5]/(dA_fit*dwTb_fit)
corr_Bphi = covm[2][3]/(dB_fit*dphi_fit)
corr_BwTa = covm[2][4]/(dB_fit*dwTa_fit)
corr_BwTb = covm[2][5]/(dB_fit*dwTb_fit)
corr_phiwTa = covm[3][4]/(dphi_fit*dwTa_fit)
corr_phiwTb = covm[3][5]/(dphi_fit*dwTb_fit)
corr_wTawTb = covm[4][5]/(dwTa_fit*dwTb_fit)
corm = np.zeros((6,6))
for i in range(6):
    for j in range (6):
        corm[i][j] = covm[i][j]/covm[i][i]

print('Matrice di correlazione:\n', corm)    
print('Covarianza normalizzata fA:', corr_fA)
print('Covarianza normalizzata fB:', corr_fB)
print('Covarianza normalizzata fphi:', corr_fphi)
print('Covarianza normalizzata fwTa:', corr_fwTa)
print('Covarianza normalizzata fwTb:', corr_fwTb)
print('Covarianza normalizzata AB:', corr_AB)
print('Covarianza normalizzata Aphi:',corr_Aphi)
print('Covarianza normalizzata AwTa:',corr_AwTa)
print('Covarianza normalizzata AwTb:',corr_AwTb)
print('Covarianza normalizzata Bphi:',corr_Bphi)
print('Covarianza normalizzata BwTa:',corr_BwTa)
print('Covarianza normalizzata BwTb:',corr_BwTb)
print('Covarianza normalizzata phiwTa:',corr_phiwTa)
print('Covarianza normalizzata phiwTb:',corr_phiwTb)
print('Covarianza normalizzata wTawTb:',corr_wTawTb)

#Plot DV vs t
if tex:
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

tt = np.linspace(min(st), max(st), 2000)
fig1,(ax1, ax2) = plt.subplots(2,1, True, gridspec_kw={'wspace':0.05,
     'hspace':0.05, 'height_ratios': [3, 1]})
ax1.set_ylabel('Differenza di potenziale $\Delta V$ [V]', y=0.6)
ax1.grid(color = 'gray', ls = '--', alpha=0.7)
ax1.errorbar(st*1e3, sDV, dDV, dt, 'ko', ms=1.5, elinewidth=1., capsize=1.5,
              ls='',label='data', zorder=5)
ax1.plot(tt*1e3, cas(tt, *pars), c='gray',
        label='fit\n $\chi^2 = 5\\times 10^6/%d$' %ndof, zorder=10, alpha=0.7)
if tick:
    ax1.yaxis.set_major_locator(plt.MultipleLocator(1.))
    ax1.yaxis.set_minor_locator(plt.MultipleLocator(0.2))
ax1.tick_params(direction='in', length=5, width=1., top=True, right=True)
ax1.tick_params(which='minor', direction='in', width=1., top=True, right=True)
legend = ax1.legend(loc ='best', framealpha = 0.3)
legend.set_zorder(100)

ax2.set_xlabel('Tempo $t$ [ms]', x=0.92)
ax2.set_ylabel('Residui')
ax2.axhline(0, c='r', zorder=10)
ax2.errorbar(st*1e3, resnorm, None, None, 'ko', elinewidth = 0.7,
              capsize=0.7, ms=1., ls='--', lw=0.8, zorder=0)
ax2.grid(color ='gray', ls = '--', alpha=0.7)
ax2.ticklabel_format(axis='both', style='sci', scilimits=None,
                     useMathText=True)
if tick:
    ax2.xaxis.set_major_locator(plt.MultipleLocator(5))
    ax2.xaxis.set_minor_locator(plt.MultipleLocator(1))
    ax2.yaxis.set_major_locator(plt.MultipleLocator(20))
    ax2.yaxis.set_minor_locator(plt.MultipleLocator(5))
ax2.tick_params(direction='in', length=5, width=1., top=True, right=True)
ax2.tick_params(which='minor', direction='in', width=1., top=True, right=True)

# Fit con rimozione degli outliers
TT=np.array([])
VV=np.array([])
dTT=np.array([])
dVV=np.array([])

outT=np.array([])
outV=np.array([])
doutT=np.array([])
doutV=np.array([])
'Contatori in e out liers'
j=0
k=0
#tengo solo i dati che si discostano dal modello per meno
# di 3.5 deviazioni standard
soglia = 3.5
for i in range (len(sDV)):
    if (np.abs(sDV[i] - cas(st, *pars)[i])< soglia*dDV[i]): 
        TT=np.insert(TT, j, st[i])
        dTT=np.insert(dTT, j, dt[i])
        VV=np.insert(VV, j, sDV[i])
        dVV=np.insert(dVV, j, dDV[i])
        j+=1
    else:
        outT=np.insert(outT, k, st[i])
        doutT=np.insert(doutT, k, dt[i])
        outV=np.insert(outV, k, sDV[i])
        doutV=np.insert(doutV, k, dDV[i])
        k+=1

pars, covm = curve_fit(cas, TT, VV, init, dVV, absolute_sigma = False)
f_fit, A_fit, B_fit, phi_fit, wTa_fit, wTb_fit = pars
df_fit, dA_fit, dB_fit, dphi_fit, dwTa_fit, dwTb_fit = np.sqrt(covm.diagonal())
print('Parametri del fit:\n', pars)
print('Matrice di Covarianza:\n', covm)
print('f = %f +- %f Hz' % (f_fit, df_fit))
print('A = %f +- %f mV' % (A_fit, dA_fit))
print('B = %f +- %f mV' % (B_fit, dB_fit))
print('phi = %f +- %f rad' % (phi_fit, dphi_fit))
print('wTa = %f +- %f rad/s' % (wTa_fit, dwTa_fit))
print('wTb = %f +- %f rad/s' % (wTb_fit, dwTb_fit))
#Test Chi quadro per cas
ndof= len(TT) - len(pars)
normin = (VV-cas(TT, *pars))/dVV
normout = (outV-cas(outT, *pars))/doutV
chisqin = (normin**2).sum()
chirid = chisqin/ndof
sigmain = (chisqin - ndof)/np.sqrt(2*ndof)
print('Chi quadro/ndof = %f/%d [%+.1f]' % (chisqin, ndof, sigmain))
print('Chi quadro ridotto:', chirid)

# Covarianza tra i parametri di cas
corr_fA = covm[0][1]/(df_fit*dA_fit)
corr_fB = covm[0][2]/(df_fit*dB_fit)
corr_fphi = covm[0][3]/(df_fit*dphi_fit)
corr_fwTa = covm[0][4]/(df_fit*dwTa_fit)
corr_fwTb = covm[0][5]/(df_fit*dwTb_fit)
corr_AB = covm[1][2]/(dA_fit*dB_fit)
corr_Aphi = covm[1][3]/(dA_fit*dphi_fit)
corr_AwTa = covm[1][4]/(dA_fit*dwTa_fit)
corr_AwTb = covm[1][5]/(dA_fit*dwTb_fit)
corr_Bphi = covm[2][3]/(dB_fit*dphi_fit)
corr_BwTa = covm[2][4]/(dB_fit*dwTa_fit)
corr_BwTb = covm[2][5]/(dB_fit*dwTb_fit)
corr_phiwTa = covm[3][4]/(dphi_fit*dwTa_fit)
corr_phiwTb = covm[3][5]/(dphi_fit*dwTb_fit)
corr_wTawTb = covm[4][5]/(dwTa_fit*dwTb_fit)
corm = np.zeros((6,6))
for i in range(6):
    for j in range (6):
        corm[i][j] = covm[i][j]/covm[i][i]

print('Matrice di correlazione:\n', corm)    
print('Covarianza normalizzata fA:', corr_fA)
print('Covarianza normalizzata fB:', corr_fB)
print('Covarianza normalizzata fphi:', corr_fphi)
print('Covarianza normalizzata fwTa:', corr_fwTa)
print('Covarianza normalizzata fwTb:', corr_fwTb)
print('Covarianza normalizzata AB:', corr_AB)
print('Covarianza normalizzata Aphi:',corr_Aphi)
print('Covarianza normalizzata AwTa:',corr_AwTa)
print('Covarianza normalizzata AwTb:',corr_AwTb)
print('Covarianza normalizzata Bphi:',corr_Bphi)
print('Covarianza normalizzata BwTa:',corr_BwTa)
print('Covarianza normalizzata BwTa:',corr_BwTb)
print('Covarianza normalizzata phiwTa:',corr_phiwTa)
print('Covarianza normalizzata phiwTb:',corr_phiwTb)
print('Covarianza normalizzata wTawTb:',corr_wTawTb)

# Plot DV vs t con outliers di cas
fig2,(ax1, ax2) = plt.subplots(2,1, True, gridspec_kw={'wspace':0.05,
     'hspace':0.05, 'height_ratios': [3, 1]})
if(log):
    ax1.set_yscale('log')
    ax1.minorticks_on()
ax1.errorbar(TT*1e3, VV, dVV, dTT, 'ko',  ms=1.5, elinewidth=1.,
             capsize=1.5, ls='', label='data')
ax1.errorbar(outT*1e3, outV, doutV, doutT, 'gx',  ms=3, elinewidth=1.,
             capsize=1.5, ls='', label='outliers')
ax1.set_ylabel('Differenza di potenziale $\Delta V$ [V]', y=0.6)
ax1.grid(which = 'major', color = 'gray', ls = '--', alpha=0.7)
ax1.plot(tt*1e3, cas(tt, *pars), c='gray',
         label='fit \n $\chi^2 = %.f/%d$' %(chisqin, ndof), zorder=10, alpha=0.7)
if tick:
    ax1.yaxis.set_major_locator(plt.MultipleLocator(1.))
    ax1.yaxis.set_minor_locator(plt.MultipleLocator(0.2))
ax1.tick_params(direction='in', length=5, width=1., top=True, right=True)
ax1.tick_params(which='minor', direction='in', width=1., top=True, right=True)
legend = ax1.legend(loc ='best', framealpha = 0.3)
legend.set_zorder(100)

ax2.set_xlabel('Tempo $t$ [ms]', x=0.92)
ax2.set_ylabel('Residui')
ax2.axhline(0, c='r', alpha=0.7, zorder=10)
ax2.errorbar(TT*1e3, normin, None, None, 'ko', elinewidth = 0.5, capsize=0.5,
             ms=1., ls='--', lw=1., zorder=5)
ax2.errorbar(outT*1e3, normout, None, None, 'gx', elinewidth = 0.7, capsize=0.7,
             ms=3., ls='', zorder=5)
ax2.grid(color ='gray', ls = '--', alpha=0.7)
ax2.ticklabel_format(axis='both', style='sci', scilimits=None,
                     useMathText=True)
if tick:
    ax2.xaxis.set_major_locator(plt.MultipleLocator(5))
    ax2.xaxis.set_minor_locator(plt.MultipleLocator(1))
    ax2.yaxis.set_major_locator(plt.MultipleLocator(2))
    ax2.yaxis.set_minor_locator(plt.MultipleLocator(0.5))
    ax2.set_ylim(min(normin)-np.std(normin), max(normin)+np.std(normin))
ax2.tick_params(direction='in', length=5, width=1., top=True, right=True)
ax2.tick_params(which='minor', direction='in', width=1., top=True, right=True)
plt.show()