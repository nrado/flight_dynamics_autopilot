#!/usr/bin/env python

import matplotlib
matplotlib.use ('TKAgg')

from pylab import *

from scipy.integrate import odeint
from scipy.signal import lti, impulse, step


def DrawDiag_GtkAgg (X, Y, FigureNum=None, NumSubplot=None, FigTitle='No title', FigSuptitle='No suptitle', Xlabel=None, Ylabel=None, type='normal', geometry=(8, 5), render=False):
    figure (figsize=geometry, facecolor='w', edgecolor='w')
    clf ()
    suptitle (FigTitle, fontsize=12)
    suptitle ('\n\n'+FigSuptitle+'\n', fontsize=9)
    subplots_adjust (left=.16, bottom=.11, right=.96, top=.9, hspace=.2)

    for numSub in range (NumSubplot):
        subplot (NumSubplot, 1, numSub+1)
        if len(Xlabel) == NumSubplot and len(Ylabel) == NumSubplot:
            xlabel (Xlabel [numSub], fontsize=11)#, style='italic')
            ylabel (Ylabel [numSub])
        else:
            print 'Error: not enough labels for subplots'
        
        if type == 'normal':
            grid ()
            plot (X[numSub], Y[numSub])
        elif type == 'semilogx':
            grid (which='major')
            grid (which='minor')
            semilogx (X[numSub], Y[numSub])
        elif type == 'step':
            grid ()
            plot (step(Y[numSub])[0], step(Y[numSub])[1])
        elif type == 'impulse':
            grid ()
            plot (impulse (Y[numSub])[0], impulse (Y[numSub])[1])
        else:
            print 'Error: Specify plot type'
    if render:
        show ()
