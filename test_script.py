##NOTE: in order to display the histogram, the script must be run with the -i option
import numpy as np
import matplotlib
#matplotlib.use("macOSX")
#matplotlib.use("TKAgg")
import matplotlib.pyplot as plt
import os
import ROOT
import sys
from ROOT import TH1F, TH1D, TH2
#from ROOT import TCanvas
#from ROOT import gROOT
#from scipy.interpolate import CubicSpline
#from scipy.interpolate import spline
#from scipy.interpolate import InterpolatedUnivariateSpline as spline
#from matplotlib.lines import Line2D
#from array import array

# The PyEBL python class contains methods useful for handling EBL SEDs.
#import PyEBL
# The PyBlazar python class contains methods useful for handling blazar spectra
# and their associated optical depth calculations, etc.
#import PyBlazar
#from ebltable.ebl_from_model import EBL

xlist = np.linspace(1, 3, 3)
#xlist = np.linspace(-3.0, 3.0, 100)
ylist = np.linspace(1, 3.0, 3)
#ylist = np.linspace(-3.0, 3.0, 100)
X, Y = np.meshgrid(xlist, ylist)
Z = np.sqrt(X**2 + Y**2)
fig,ax=plt.subplots(1,1)
cp = ax.contourf(X, Y, Z)
fig.colorbar(cp) # Add a colorbar to a plot
ax.set_title('Filled Contours Plot')
#ax.set_xlabel('x (cm)')
ax.set_ylabel('y (cm)')
plt.show()

print Z
