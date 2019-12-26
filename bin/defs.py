# -*- coding: utf-8 -*-
"""
Created on Sat Dec 14 00:37:16 2019

@author: berna
"""
import numpy as np
from matplotlib import pyplot as plt
from collections import namedtuple

# Definizioni variabili per la serie di Fourier
t0 = -2.         # tempi iniziale
tN = 2.         # e finale
n_points = 10000 # punti equispaziati
sums = 10000     # termini nella serie

wT = 10 # frequenza angolare di taglio

# Array di appoggio per lo sviluppo in serie
tt=np.linspace(t0, tN, n_points)
c = np.zeros(len(tt))
w = np.zeros_like(c)
A = np.zeros_like(c)
Df = np.zeros_like(c)
# Tuple di appoggio per le funzioni
args=namedtuple('pars', 'f A_pp')
funcs = []
# Nota: di default frequenza =1 Hz e Ampiezza(p-p) =1 au 
def sqw(t, f=1, A_pp=1):
    fk=0
    for k in range(1, sums, 2):
        c[k] = 2./(k*np.pi)
        w[k] = k*2*np.pi*f
        fk += c[k]*np.sin(w[k]*t)
    return fk*A_pp

def trg(t, f=1, A_pp=1):
    fk=0
    for k in range(1, sums, 2):
        c[k] = (2./(k*np.pi))**2
        w[k] = k*2*np.pi*f
        fk += c[k]*np.cos(w[k]*t)
    return fk*A_pp

def fin(t, f=1, A_pp=1):
    fk = 0
    for k in range(1, sums, 2):
        c[k] = 2./(k*np.pi)
        w[k] = k*2*np.pi*f
        A[k] = 1./np.sqrt(1 + (w[k]/wT)**2)
        Df[k] = np.arctan(-w[k]/wT)
        fk += A[k]*c[k]*np.sin(w[k]*t + Df[k])
    return fk*A_pp

def plot(funcs, x, args, fit=False, log=False, tick=False, tex=False):
    # Grafico segnale su tempo
    fig, ax = plt.subplots()
    if log:
        x = np.logspace(np.log10(min(x)), np.log10(min(x)), n_points)
        ax.set_yscale('log')
        ax.set_xscale('log')
    for f in funcs:
        ax.plot(x, f(x, *args))
    ax.grid(color = 'gray', ls = '--', alpha=0.7)
    ax.set_xlabel('Tempo $t$ [s]', x=0.94)
    ax.set_ylabel('Ampiezza [arb. un]')
    ax.minorticks_on()
    if tick:
        ax.xaxis.set_major_locator(plt.MultipleLocator((tN-t0)/8.))
        ax.xaxis.set_minor_locator(plt.MultipleLocator((tN-t0)/40.))
        # ax.yaxis.set_major_locator(plt.MultipleLocator(0.1))
        # ax.yaxis.set_minor_locator(plt.MultipleLocator(0.02))
    ax.tick_params(direction='in', length=5, width=1., top=True, right=True)
    ax.tick_params(which='minor', direction='in', width=1., top=True, right=True)
    # Impostazioni opzionali per grafici in LaTeX
    if tex:
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
    plt.show()