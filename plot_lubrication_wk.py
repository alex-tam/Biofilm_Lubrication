#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot numerical solutions to the 1D extensional flow model

Created on Tue Jul 3 14:20:25 2018
@author: Alex Tam
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib2tikz as m2t

# %% Import numerical data
times = np.arange(1, 13002, 1000) # time steps at which to import solutions
r = np.genfromtxt('r.csv', delimiter = ',') # import numerical r

# %% Plot numerical solution
plt.rc('text', usetex=True) # allow LaTeX in labels
plt.rc('xtick', labelsize=20); plt.rc('ytick', labelsize=20) # set tick label size

for i in range(0, np.size(times)):
    h = np.genfromtxt('biofilm_height-'+str(times[i])+'.csv', delimiter = ',')
    if (i == 0):
        plt.plot(r, h, 'k--'); plt.axis([0, 10, 0, 1.5]) # set (x,y) axis limits
    else:
        plt.plot(r, h, 'k'); plt.axis([0, 10, 0, 1.5]) # set (x,y) axis limits
plt.annotate("", xy=(9, 1), xytext=(0.25, 0.25),arrowprops=dict(arrowstyle="->")) # manually draw arrow
plt.text(9, 1, r'\(t\)', fontsize=20)
plt.xlabel(r'\(r\)', fontsize=20); plt.ylabel(r'\(h(r,t)\)', fontsize=20); plt.tight_layout(); 
m2t.save('height_lubrication_wk.tex')
plt.show(); plt.figure()

for i in range(0, np.size(times)):
    J = np.genfromtxt('source-'+str(times[i])+'.csv', delimiter = ',')
    if (i == 0):
        plt.plot(r, J, 'k--'); plt.axis([0, 10, -0.1, 0.15]) # set (x,y) axis limits
    else:
        plt.plot(r, J, 'k'); plt.axis([0, 10, -0.1, 0.15]) # set (x,y) axis limits
plt.annotate("", xy=(9, -0.04), xytext=(0.25, 0.075),arrowprops=dict(arrowstyle="->")) # manually draw arrow
plt.text(9, -0.04, r'\(t\)', fontsize=20)
plt.xlabel(r'\(r\)', fontsize=20); plt.ylabel(r'\(J(r,t)\)', fontsize=20); plt.tight_layout();
m2t.save('source_lubrication_wk.tex')
plt.show()