#!/usr/bin/env python

import numpy as np
from scipy import interpolate
from scipy import integrate
import math
#import pygsl.monte
#import pygsl.rng
# import pygsl._numobj as Numeric
import types
import numbers
import linecache

import matplotlib.pyplot as plt
#import matplotlib.ticker as ticker

#from matplotlib import rc
#rc('font',size=16)
#rc('xtick.major',size=7)
#rc('xtick.minor',size=4)
#rc('ytick.major',size=7)
#rc('ytick.minor',size=4)

# Thomson cross-section in units of meters-squared
thomson_cross_sec = 6.65e-25/1.e4
# Speed of light in units of meters per second
speed_of_light = 2.998e8
# Planck's constant in units of eV-seconds
planck_const = 4.136e-15
# Electron mass in units of kilograms
electron_mass = 9.109e-31
# eV per Joule
eVPerJoule = 1./1.6e-19
# Meters per Mpc
meterPerMpc = 3.08568025e22

# Hubble's constant in units of meters per second per Mpc
hubble_const = 70.e3
# Hubble radius in units of Mpc
hubble_radius = speed_of_light/hubble_const
# Cosmological matter density
omega_matter = 0.27
# Cosmological radiation density
omega_rad = 0.0
# Cosmological constant density
omega_lambda = 0.73

class ebl:
    
    def __init__(self, fileName, linenum):
        
        self.eblSEDArr = []
        self.eblTau = []
        self.enGamma = []
        self.redshift = 0.1
        
        # Hubble's constant in units of meters per second per Mpc
        self.hubble_const = 70.e3
        # Hubble radius in units of Mpc
        self.hubble_radius = speed_of_light/self.hubble_const
        # Cosmological matter density
        self.omega_matter = 0.27
        # Cosmological radiation density
        self.omega_radiation = 0.0
        # Cosmological constant density
        self.omega_lambda = 0.73
        
        self.eblSED = self.loadEBLSED(fileName, linenum)
    
    ############################################################################
    #                             Load an EBL SED                              #
    ############################################################################
    def loadEBLSED(self, fileName, linenum):
        
        FILEIN = open(fileName,"r")
        
        self.eblSEDArr = []
        
        line = linecache.getline(fileName, linenum)
        self.eblSEDArr = line.split()
        #print self.eblSEDArr

        #hard-coded the knot positions in lambda to make things easy
        #lamb=[10]
        #lamb=[0.178,  0.316,  0.563,  1,  1.78,  3.16,  5.62,  10,  17.8,  31.6, 56.3, 100]
        #lamb=[5,  6.05764,  7.339,  8.8914,  10.7722,  13.0508,  15.8114,  19.1559,  23.2079,  28.1171,  34.0646,  41.2702]
        #lamb=[1,  1.77828,  3.16228,  5.62341,  10,  17.7828,  31.6228,  56.2341, 100,  177.828,  316.228,  562.341]
        #lamb=[1,  84.25,  167.5,  250.75,  334,  417.25,  500.5,  583.75,  667,  750.25,  833.5,  916.75] 
        #lamb=[1,  1.77828,  3.16228,  5.62341,  10,  17.7828,  31.6228,  56.2341, 100,  177.828,  316.228,  562.341]
        lamb= np.logspace(np.log10(1),np.log10(100),12, endpoint=True)
        return interpolate.UnivariateSpline(lamb, self.eblSEDArr, k=2, s=0)
        
    ############################################################################
    #          ---------------------------------------------------             #
    #          | Methods for gamma-ray optical depth calculation |             #
    #          ---------------------------------------------------             #
    ############################################################################
        
    ############################################################################
    # EBL photon energy threshold for pair production with a gamma-ray with    #
    # energy enGamma.  Units are eV.                                           #
    ############################################################################
    def enThresh(self, enGamma, mu):
        return (2.*(electron_mass*speed_of_light**2)**2. * eVPerJoule**2) / (enGamma*(1.-mu))
    
    def beta(self, enGamma,enEBL, mu):
        if self.enThresh(enGamma,mu) < enEBL:
            return np.sqrt(1. - self.enThresh(enGamma,mu)/enEBL)
        else:
            return 0.
    
    ############################################################################
    # Thomson cross-section for energies enGamma and enEBL with the angle      #
    # between photons being mu=cos(theta)                                      #
    ############################################################################
    def cross_sec(self, enGamma, enEBL, mu):
        prefactor = (3.*thomson_cross_sec/16.) * (1. - self.beta(enGamma,enEBL,mu)**2.)
        factor1 = 2.*self.beta(enGamma,enEBL,mu)*(self.beta(enGamma,enEBL,mu)**2 - 2.)
        factor2 = (3. - self.beta(enGamma,enEBL,mu)**4)*np.log((1.+self.beta(enGamma,enEBL,mu))/(1.-self.beta(enGamma,enEBL,mu)))
        return prefactor*(factor1 + factor2)
    
    def cosmo_line_element(self, z):
        factor = np.sqrt( (1.+z)**2*(omega_matter*z+1.) + z*(2.+z)*((1.+z)**2*omega_rad - omega_lambda) )
        return hubble_radius/((1.+z)*factor)
    
    ############################################################################
    #                       EBL photon number density                          #
    ############################################################################
    def ebl_number_density(self, enEBL):
        prefactor = 4.*math.pi/(speed_of_light*enEBL**2)*10.**-9*eVPerJoule
        return prefactor*self.eblSED(planck_const*speed_of_light*10.**6/enEBL)
    
    ############################################################################
    #            Integrand for gamma-ray optical depth calculation             #
    ############################################################################
    # Note: the ordering of the arguments here is very important.  The variable which is integrated
    # over first must be the first argument, the second variable integrated over must be the second
    # argument, and so on.  Parameters that are held fixed, and get entered into the integreal as
    # 'args' must come after all of the integration variables.
    def tau_integrand(self, enEBL, mu, z, enGamma, fEvo):
#        print "z: %f  mu: %f  enEBL: %g  enThresh: %g  num_dens: %g  cross_sec: %g" % (z, mu, enEBL, self.enThresh(enGamma,mu), self.ebl_number_density(enEBL), self.cross_sec(enGamma, enEBL, mu)) 
#        print "enGamma: %f  fEvo: %f" % (enGamma, fEvo)
        return ( meterPerMpc*self.cosmo_line_element(z)*((1.-mu)/2.)*(1.+z)*(1.+z)**(3.-fEvo)
                *self.ebl_number_density(enEBL)*self.cross_sec(enGamma*(1+z), enEBL, mu) )
    
    # The monte carlo integration method makes use of PyGSL, which makes use of the GNU
    # Scientific Libraries (GSL).  For this reason, when the monte carlo integration technique
    # is used, the arguments to the integrand have to be handled differently.  Rather than
    # specifying each integration parameter and fixed paramter separately, two arrays are
    # given.  The first array contains the integration variables (no values are actaully passsed)
    # and the second array contains fixed paramters.
    def tau_integrand_monte(self, integrand_vars, params):
        return ( meterPerMpc*self.cosmo_line_element(integrand_vars[2])*((1.-integrand_vars[1])/2.)
                *(1.+integrand_vars[2])*(1.+integrand_vars[2])**(3.-params[1])
                *self.ebl_number_density(integrand_vars[0])*self.cross_sec(params[0]*(1.+integrand_vars[2]), integrand_vars[0], integrand_vars[1]) )

    ############################################################################
    #         Calculate the gamma-ray opacity for the loaded EBL SED           #
    ############################################################################
    def calcGammaAbsorp(self, z, enGamma=1.e12, fEvo=0.,method='tplquad'):
        
        if isinstance(enGamma, types.ListType) == False:
            enGamma = [enGamma]
        
        tauResults = np.zeros((len(enGamma),2))
        
        for enIndex in range(len(enGamma)):        
            # Numerical integration using scipy.
            if method == 'tplquad':
                tau = integrate.tplquad(self.tau_integrand, 0., z, lambda mu: -1., lambda mu: 1.,
                                       lambda enEBL, mu: self.enThresh((1.+z)*enGamma[enIndex],mu),
                                       lambda enEBL, mu: (planck_const*speed_of_light/(0.1*10.**-6)),
                                       args=(enGamma[enIndex], fEvo), epsabs=0.1, epsrel=0.1)
            
            # Monte carlo integration using GSL through PyGSL.
            elif method == 'montecarlo':
                params = [enGamma[enIndex], fEvo]
                tau_integrand_gsl = pygsl.monte.gsl_monte_function(self.tau_integrand_monte, params, 3)
                
                r = pygsl.rng.mt19937_1999()
                
                # low_lims = [0., -1., self.enThresh((1.+z)*enGamma,mu)]
                low_lims = [self.enThresh((1.+z)*enGamma[enIndex],-1), -1., 0.]
                up_lims = [(planck_const*speed_of_light/(0.1*10.**-6)), 1., z]
                
                #print low_lims
                #print up_lims
                
                s = pygsl.monte.miser(3)
                s.init()
                
                tau = s.integrate(tau_integrand_gsl, low_lims, up_lims, 6000, r)

                #if tau[0] > 45:
                #    break

            else:
                print("ERROR: invalid integration method specified.  Must be 'tplquad' or 'montecarlo' ")
                tau = [0,0]
                break
        
            #print"Redshift: %.3f;   Energy: %.2f TeV;   Opacity: %.4f" % (
            #    z, enGamma[enIndex]/10.**12, tau[0] )
            
            tauResults[enIndex] = tau
        
        #print("Done!\n")
        
        return tauResults
        
    
    ############################################################################
    #      ---------------------------------------------------------           #
    #      | Methods for plotting the EBL SED, optical depth, etc. |           #
    #      ---------------------------------------------------------           #
    ############################################################################    
    
    ############################################################################
    #                        Plot the loaded EBL SED                           #
    ############################################################################
    def plotEBLSED(self):
        fig = plt.figure()
        # ax = plt.axes([0.15,0.14,0.82,0.825])
        ax = plt.axes([0.165,0.14,0.8,0.82])
        
        log10x = np.linspace(-0.75,1.5,1000)
        x = 10.**log10x

        # x = np.linspace(-1.,1.,1000)
        y = self.eblSED(x)
        
        # y = (planck_const*speed_of_light/x)*self.ebl_number_density(planck_const*speed_of_light*10.**6/x)
        # y = self.cross_sec(1.e12, planck_const*speed_of_light*10.**6/x, 0.)
        # y = []
        # for i in range(len(x)):
        #     # y.append(self.cross_sec(1.e12, planck_const*speed_of_light*10.**6/x[i], -0.5)/thomson_cross_sec)
        #     y.append(self.tau_integrand(0.05,x[i],planck_const*speed_of_light*10.**6/1.,1.e12,0.))
            
            # print "%f %g" % (x[i],y[i])
        
        plt.loglog(x,y,linewidth=2.)

        #lower limits from Dwek&Krennrich, loosened by 2 sigma
        limxl=[0.153, 0.1595, 0.231, 0.2365, 0.3, 0.45, 0.61, 0.81, 1.25, 1.25, 1.6, 2.12, 2.2, 3.6, 4.5, 15, 16, 24]
        limyl=[1.03-2*0.15, 3.75-2*1.25, 2.25-2*0.32, 3.6-2*0.5, 3.7-2*0.7, 6.1-2*1.8, 7.4-2*1.5, 9.3-2*1.6, 11.5-2*1.3, 11.7-2*2.6, 11.5-2*1.5, 10.0-2*0.8, 9.0-2*1.2, 5.6-2*1.0, 4.4-2*0.8, 1.9-2*0.4, 2.2-2*0.2, 2.86-2*0.19]  
        plt.loglog(limxl, limyl, 'ro')

        #upper limits from Dwek&Krennrich, loosened by 2 sigma
        limxu=[0.3, 0.4, 0.44, 0.5115, 0.55, 0.64, 0.814, 1.25, 2.2, 3.5, 4.9]     
        limyu=[18+2*12, 26+2*10, 7.9+2*4.0, 30+2*9, 55+2*27, 7.7+2*5.8, 57+2*32, 21+2*15, 20+2*6, 13.3+2*2.8, 22+2*12] 
        plt.loglog(limxu, limyu, 'ro')

        plt.xlabel("$\lambda \, [\mu$m]",fontsize=22)
        # plt.ylabel(r"$\nu I_{\nu} \, [\mathrm{nW} \mathrm{m}^{-2} \mathrm{sr}^{-1}]$",fontsize=18)
        # plt.ylabel("$\\nu I_{\\nu} \,$[nW m$^{-2}$sr$^{-1}$]",fontsize=18)
        plt.ylabel("$\\nu I_{\\nu}$ [nW$\,$m$^{-2}$sr$^{-1}$]",fontsize=18)
        
        for label in ax.xaxis.get_ticklabels():
            label.set_fontsize(20)
        for label in ax.yaxis.get_ticklabels():
            label.set_fontsize(20)
        
        #ax.set_xlim(0.1,1000)
        #ax.set_ylim(1.,20.)
        
        plt.draw()
        plt.show()
#        plt.savefig("EBLmodel.png")
        return fig
    
    ############################################################################
    #                   Plot the EBL photon energy density                     #
    ############################################################################
    def plotEBLNumDens(self):
        fig = plt.figure()
        ax = plt.axes([0.165,0.14,0.8,0.82])
        
        #log10x = np.linspace(-1.,3.,1000)
        log10x = np.linspace(0.5,2,1000)
        x = 10.**log10x
        y = (planck_const*speed_of_light/x)*self.ebl_number_density(planck_const*speed_of_light*10.**6/x)
        
        plt.loglog(x,y,linewidth=2.)
        
        plt.xlabel("$\lambda \, [\mu$m]",fontsize=22)
        plt.ylabel(r"$\epsilon n_{\epsilon}\,$[cm$^{-3}$]",fontsize=22)
        
        for label in ax.xaxis.get_ticklabels():
            label.set_fontsize(20)
        for label in ax.yaxis.get_ticklabels():
            label.set_fontsize(20)
        
        #ax.set_xlim(0.1,1000)
        #ax.set_ylim(1.,20.)
        
        plt.draw()
        plt.show()

    ############################################################################
    #                  Plot the photon-photon cross-section                    #
    ############################################################################
    def plotCrossSec(self, enGamma):
        fig = plt.figure()
        ax = plt.axes([0.15,0.14,0.82,0.825])
        
        log10x = np.linspace(-1.,3.,1000)
        x = 10.**log10x
        y = []
        for i in range(len(x)):
            y.append(self.cross_sec(enGamma, planck_const*speed_of_light*10.**6/x[i], -0.5)/thomson_cross_sec)
        
        plt.loglog(x,y,linewidth=2.)
        
        plt.xlabel("$\lambda \, [\mu$m]",fontsize=22)
        plt.ylabel(r"$\sigma_{\gamma\gamma}$ [$\sigma_T$]",fontsize=22)
        
        for label in ax.xaxis.get_ticklabels():
            label.set_fontsize(20)
        for label in ax.yaxis.get_ticklabels():
            label.set_fontsize(20)
        
        #ax.set_xlim(0.1,1000)
        #ax.set_ylim(1.,20.)
        
        plt.draw()
        plt.show()

    ############################################################################
    #             Plot the integrand of the optical depth integral             #
    ############################################################################
    def plotTauIntegrand(self, z=0.1, eblWavelength=1., gammaEn=1.e12, mu=0., fEvo=0.):
        fig = plt.figure()
        ax = plt.axes([0.15,0.14,0.82,0.825])
        
        if isinstance(z, numbers.Number) == True:
            z = [z]
        if isinstance(eblWavelength, numbers.Number) == True:
            eblWavelength = [eblWavelength]
        if isinstance(enGamma, numbers.Number) == True:
            enGamma = [enGamma]
        if isinstance(mu, numbers.Number) == True:
            mu = [mu]
        if isinstance(fEvo, numbers.Number) == True:
            fEvo = [fEvo]
        
        # This if state performs an exclusive OR (XOR) to verify that only one of the
        # input variables is an array.  The XOR operator (^) is a bitwise operator but
        # also works on booleans and integers.  For this reason, the results of the array
        # length comparisons are cast as booleans.
        if ( bool(len(z) > 1) ^ bool(len(eblWavelength) > 1) ^ bool(len(enGamma) > 1) ^
             bool(len(mu) > 1) ^ bool(len(fEvo) > 1) == False ):
            print("ERROR: One, and only one, of the input values may be an array with a \
                   length greater than 1.  The other input values should be scalars.")
            return
        
        y = []
        for zTmp in z:
            for eblWavelengthTmp in eblWavelength:
                for enGammaTmp in enGamma:
                    for muTmp in mu:
                        for fEvoTmp in fEvo:
                            y.append(self.tau_integrand(planck_const*speed_of_light*10.**6/eblWavelengthTmp,
                                                        muTmp,zTmp,gammaEnTmp,fEvoTmp))
            
            # print "%f %g" % (x[i],y[i])
        
        plt.loglog(x,y,linewidth=2.)
        
        plt.xlabel("$\lambda \, [\mu$m]",fontsize=22)
        plt.ylabel(r"$\nu I_{\nu} \, [\mathrm{nW} \, \mathrm{m}^{-2} \, \mathrm{sr}^{-1}]$",fontsize=22)
        
        for label in ax.xaxis.get_ticklabels():
            label.set_fontsize(20)
        for label in ax.yaxis.get_ticklabels():
            label.set_fontsize(20)
        
        #ax.set_xlim(0.1,1000)
        #ax.set_ylim(1.,20.)
        
        plt.draw()
        plt.show()
