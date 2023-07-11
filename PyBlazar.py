#!/usr/bin/env python

import cPickle
import numpy as np
import matplotlib.pyplot as plt
#from matplotlib import rc
#rc('font',size=18)
#rc('xtick.major',size=7)
#rc('xtick.minor',size=4)
#rc('ytick.major',size=7)
#rc('ytick.minor',size=4)

import PyEBL

class blazar:
    
    def __init__(self):
        
        self.eblSED = []
        self.eblTau = []
        self.eblTauErr = []
        self.energies = []
        self.energiesTeV = []
        self.specTeV = []
        self.redshift = 0.1
        self.intrSpecTeV = []
    
    ############################################################################
    #                             Load an EBL SED                              #
    ############################################################################
    def loadEBLSED(self, fileName, linenum):

        self.eblModel = PyEBL.ebl(fileName, linenum)
        
        #FILEIN = open(fileName,"r")
        #
        #i = 0
        #for line in FILEIN:
        #    lambdaTmp, eblIntensTmp = line.split()
        #    
        #    self.eblSED.append([])
        #    self.eblSED[i].append(float(lambdaTmp))
        #    self.eblSED[i].append(float(eblIntensTmp))
        #    
        #    i += 1
    
    ############################################################################
    #                         Load a blazar spectrum                           #
    ############################################################################
    def loadTeVSpec(self, fileName='', specTeV=[]):
        
        self.energies = []
        self.specTeV = []
        
        if fileName == '' and np.shape(specTeV)[0] > 0 and np.shape(specTeV)[1] == 3:
            self.specTeV = specTeV
            for index in range(len(specTeV)):
                self.energies.append(specTeV[index][0]*10**12)
                
        elif fileName != '':
            
            FILEIN = open(fileName,"r")
            
            i = 0
            for line in FILEIN:
                enTmp, fluxTmp, fluxTmpErr = line.split()
                
                self.energies.append(float(enTmp)*10**12)
                
                self.specTeV.append([])
                self.specTeV[i].append(float(enTmp))
                self.specTeV[i].append(float(fluxTmp))
                self.specTeV[i].append(float(fluxTmpErr))
                
                i += 1
        
        else:
            print("Please specify a file or pass in an array with spectral data.\n")
    
        self.energiesTeV = [energy/1.e12 for energy in self.energies]
    
        self.specTeV = np.array(self.specTeV)
        
    ############################################################################
    #            Load the gamma-ray opacity for a given EBL SED                #
    ############################################################################
    def loadGammaAbsorp(self, fileName):
        
        self.eblTau = []
        
        FILEIN = open(fileName,"r")
        
        i = 0
        for line in FILEIN:
            redshiftTmp, enTmp, tauTmp = line.split()
            
            # self.eblTau.append([])
            # self.eblTau[i].append(float(enTmp))
            # self.eblTau[i].append(float(tauTmp))
            self.eblTau.append(float(tauTmp))
            
            i += 1
    
    ############################################################################
    # Check that the energies of the loaded gamma-ray opacity match the        #
    # energies of the loaded blazar spectrum.                                  #
    ############################################################################        
    def validateTeVGammaAbsorp(self):
        
        if len(eblTau) != len(specTeV):
            print("ERROR: EBL opacity and blazar spectrum have a different number of energies!\n")
            
            return False
        else:
            for i in range(len(specTeV)):
                eps = specTeV[i][0] - eblTau[i][0]
                if eps > 1.e-2:
                    print("ERROR: energies for blazar spectrum and gamma-ray opacity do not match!\n")
                    
                    return False
        
        return True
                    
    
    ############################################################################
    #         Calculate the gamma-ray opacity for the loaded EBL SED           #
    ############################################################################
    def calcGammaAbsorp(self, fEvo=0.,method='montecarlo'):
        self.eblTau = []
        self.eblTauErr = []

        results = self.eblModel.calcGammaAbsorp(self.redshift, self.energies, fEvo=fEvo, method=method)
        
        self.eblTau = [row[0] for row in results]
        self.eblTauErr = [row[1] for row in results]
        
        self.intrSpecTeV = []
        for index in range(len(self.specTeV)):
            self.intrSpecTeV.append([])
            self.intrSpecTeV[index].append(self.specTeV[index][0])
            self.intrSpecTeV[index].append(np.exp(self.eblTau[index])*self.specTeV[index][1])
            self.intrSpecTeV[index].append(np.exp(self.eblTau[index])*self.specTeV[index][2])

        self.intrSpecTeV = np.array(self.intrSpecTeV)
    
    ############################################################################
    #        ------------------------------------------------------            #
    #        | Methods for plotting optical depth, spectrum, etc. |            #
    #        ------------------------------------------------------            #
    ############################################################################
    def plotTau(self,newFig=True):
        if newFig == True:
            fig = plt.figure()
            # ax = plt.axes([0.15,0.15,0.82,0.815])
            ax = plt.axes([0.18,0.15,0.805,0.81])
        
        energies = [energy/1.e12 for energy in self.energies]
        
        plt.loglog(energies,self.eblTau,linewidth=2.)
        
        plt.xlabel("$E_{\gamma} \,$[TeV]",fontsize=22)
        plt.ylabel(r"$\tau_{\gamma\gamma}$",fontsize=22)
        
        # plt.xlim(energies[0,0]/2.,energies[0,-1]*2.)
        plt.xlim(energies[0]/2.,energies[-1]*2.)
        
        if newFig == True:
            for label in ax.xaxis.get_ticklabels():
                label.set_fontsize(20)
            for label in ax.yaxis.get_ticklabels():
                label.set_fontsize(20)
        
        #ax.set_xlim(0.1,1000)
        #ax.set_ylim(1.,20.)
        
        plt.draw()
        plt.show()
        
    def plotSpectrum(self, ESq=False):
        self.fig = plt.figure(figsize=(7,6))
        # ax = plt.axes([0.16,0.12,0.82,0.83])
        ax = plt.axes([0.18,0.15,0.805,0.81])
        ax.set_xscale('log')
        ax.set_yscale('log')
        
        if ESq == True:
            plt.errorbar(self.specTeV[:,0], self.specTeV[:,1]*self.specTeV[:,0]**2, yerr=self.specTeV[:,2]*self.specTeV[:,0]**2,
                         linestyle='None', marker='o', markeredgewidth=1.5, color='black', markerfacecolor='green')
            plt.ylabel('$E^2 dN/dE$ [TeV$\,$cm$^{-2}$s$^{-1}$]')
        else:
            plt.errorbar(self.specTeV[:,0], self.specTeV[:,1], yerr=self.specTeV[:,2], linestyle='None', marker='o',
                         markeredgewidth=1.5, color='black', markerfacecolor='green')
            plt.ylabel('$dN/dE$ [ph$\,$cm$^{-2}$s$^{-1}$]')
        
        plt.xlabel('Energy [TeV]')
        
        plt.xlim(self.specTeV[0,0]/2., self.specTeV[-1,0]*2.)
        
        plt.draw()
        plt.show()
        
        # return fig
    
    def plotIntrSpectrum(self, ESq=False):
        self.fig = plt.figure(figsize=(7,6))
        # ax = plt.axes([0.16,0.12,0.82,0.83])
        ax = plt.axes([0.18,0.15,0.805,0.81])
        ax.set_xscale('log')
        ax.set_yscale('log')
        
        if ESq == True:
            plt.errorbar(self.intrSpecTeV[:,0], self.intrSpecTeV[:,1]*self.intrSpecTeV[:,0]**2, yerr=self.intrSpecTeV[:,2]*self.intrSpecTeV[:,0]**2,
                         linestyle='None', marker='o', markeredgewidth=1.5, color='black', markerfacecolor='green')
            plt.ylabel('$E^2 dN/dE$ [TeV$\,$cm$^{-2}$s$^{-1}$]')
        else:
            plt.errorbar(self.intrSpecTeV[:,0], self.intrSpecTeV[:,1], yerr=self.intrSpecTeV[:,2], linestyle='None', marker='o',
                         markeredgewidth=1.5, color='black', markerfacecolor='green')
            plt.ylabel('$dN/dE$ [ph$\,$cm$^{-2}$s$^{-1}$]')
        
        plt.xlabel('Energy [TeV]')
        
        plt.xlim(self.specTeV[0,0]/2., self.specTeV[-1,0]*2.)
        
        plt.draw()
        plt.show()

        # return fig
        
    
    ############################################################################
    #              -----------------------------------------------             #
    #              | Method for saving the blazar python object |              #
    #              -----------------------------------------------             #
    ############################################################################
    def save(self,fileName):
        
        pklFile = open(fileName,"wb")
        cPickle.dump(self,pklFile)
        pklFile.close()
