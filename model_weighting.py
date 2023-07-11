##NOTE: in order to display the histogram, the script must be run with the -i option
import numpy as np
import matplotlib.pyplot as plt 
import os
import ROOT
import sys
from ROOT import TH1F, TH1D, TH2
from ROOT import TCanvas, TPad, TFile, TPaveLabel, TPaveText, TFile, TLine, TGraph, TAttFill, TColor, TLegend
from ROOT import gROOT
from scipy.interpolate import spline
from matplotlib.lines import Line2D

# The PyEBL python class contains methods useful for handling EBL SEDs.
import PyEBL
# The PyBlazar python class contains methods useful for handling blazar spectra
# and their associated optical depth calculations, etc.
import PyBlazar
from ebltable.ebl_from_model import EBL 
plt.ion()

'''
hist_1 = TH1F('hist_1', 'Intensity @ #lambda = 1#mum', 50, 1, 32)
hist_2 = TH1F('hist_1', 'Intensity @ #lambda = 84#mum', 50, 1, 32)
hist_3 = TH1F('hist_1', 'Intensity @ #lambda = 500#mum', 50, 1, 32)
hist_4 = TH1F('hist_1', 'Intensity @ #lambda = 833#mum', 50, 1, 32)
'''

# FIRST!! Retrieve Model lambda-Intensities
namewithdir = '/home/mateo/Documents/HAWC/EBL/EBLmodel_likelihood/mrk421_gp_SSC_ECPL_EBL_loglike_Fixed_FromLimits/EBLintensities_knotfree_logScale_1e4sample_fromUpperLowerBounds_includingStandardModels.dat'
f_modelIntensities=open(namewithdir, 'r')


d_modelIntensities = f_modelIntensities.readlines()

#reader = csv.reader(f, delimiter=' ', skipinitialspace=True)
#first_row = next(reader)
#num_cols = len(first_row)

#lamb=[1,  1.77828,  3.16228,  5.62341,  10,  17.7828,  31.6228,  56.2341, 100,  177.828,  316.228,  562.341]
nmodels=10000
#lamb=[5,  6.05764,  7.339,  8.8914,  10.7722,  13.0508,  15.8114,  19.1559,  23.2079,  28.1171,  34.0646,  41.2702]
#lamb=[1,  84.25,  167.5,  250.75,  334,  417.25,  500.5,  583.75,  667,  750.25,  833.5,  916.75]
lamb= np.logspace(np.log10(1),np.log10(1500),12, endpoint=False)

#### Define array of hisograms for intensity values and Tgraphs for areas

#TH1F* h_intensity_fitchi2[12];

hist_array=[]
#hist = TH1F('hist', 'Intensity', 10, -1, 32)

graph_array=[]



for l in range(len(lamb)):
	hist_name="hIntensity@ #lambda={0}#mum".format(str(lamb[l]))
	#print hist_name	
	hist = TH1F(hist_name, hist_name, 40, -1, 15)
	#hist[l].append[ROOT.TH1F("myHisto","myHisto",64,-4,4)]
	hist_array.append(hist)



Intensity=np.empty([nmodels,len(lamb)])

## Create array with intensities for each model (first index) and wavelenght (second index) 
i=-1
for x in d_modelIntensities: 
	#Intensity.append([])
        i=i+1;
	for s in range(len(lamb)):
		I=float(x.split('  ')[s])
		Intensity[i][s]=I



## SECOND!! Create array with weights
f_results=open('/home/mateo/Documents/HAWC/EBL/EBLmodel_likelihood/mrk421_gp_SSC_ECPL_EBL_loglike_Fixed_FromLimits/median/Model_BIC_Categories_Ordered_ProbValues_Posterior.txt', 'r')
#f_results=open('/home/mateo/Documents/HAWC/EBL/EBLmodel_likelihood/mrk421_gp_SSC_ECPL_EBL_loglike_Fixed_AllRange_Subset_1e3/Model_BIC_Categories_Ordered_ProbValues_fromFits_gaus.txt', 'r')
model_probValues = f_results.readlines()

#Intensity_w=[]## weight for each model (1st index) and each lamba (2nd index)
Intensity_w=np.empty([nmodels,len(lamb)])


total_points = Intensity.size
Int_array = np.empty([total_points])
Int_array_spline = np.empty([nmodels])
lambda_array = np.empty([total_points])
#sigma_array = np.empty([total_points,total_points])
#cate_array = np.empty([total_points])
cate_array = np.empty([nmodels])
k=-1
i=-1
for x in model_probValues:
	model= x.split('  ')[0]	
	cate= float(x.split('  ')[1])
	i=i+1
	cate_array[i]=cate
	#if(cate>0.01):
	#	print model," ",cate
	for j in range(len(lamb)):
		k=k+1
                Int_array[k]=Intensity[i][j]
                lambda_array[k]=lamb[j]
		#Intensity_w[i][j]=np.exp(-cate)
		Intensity_w[i][j]=cate

'''
for x in d_results:
	sigma= float(x.split(' ')[1])
	print sigma
	for i in range(total_points):
		for j in range(total_points):
		#	print i," ",j," ",sigma
			sigma_array[i][j]=sigma
'''



#plot intensity weights for each lamba:


#hist_5.GetXAxis().SetTitle("bla")
#for i in range(len(lamb)):


lower_bound_16=np.zeros(shape=(len(lamb),hist_array[j].GetSize()))
lower_containment_16=np.zeros(shape=(len(lamb),hist_array[j].GetSize()))

lower_bound_05=np.zeros(shape=(len(lamb),hist_array[j].GetSize()))
lower_containment_05=np.zeros(shape=(len(lamb),hist_array[j].GetSize()))


higher_bound_16=np.zeros(shape=(len(lamb),hist_array[j].GetSize()))
higher_containment_16=np.zeros(shape=(len(lamb),hist_array[j].GetSize()))

higher_bound_05=np.zeros(shape=(len(lamb),hist_array[j].GetSize()))
higher_containment_05=np.zeros(shape=(len(lamb),hist_array[j].GetSize()))


low_bound_16=[]
high_bound_16=[]

low_bound_05=[]
high_bound_05=[]

line_low_array=[]
line_high_array=[]
for j in range(len(lamb)):
	I=[]
	
	for i in range(len(model_probValues)):

		I.append(Intensity[i][j])
		hist_array[j].Fill(Intensity[i][j],Intensity_w[i][j])
	#First, lower bounds
	for l in range(hist_array[j].GetSize()):

		ratio =  hist_array[j].Integral(0,l)/hist_array[j].Integral()
		if (ratio > 0.16):
			#print hist_array[j].GetBinCenter(l-1)
			lower_bound_16[j][l] = hist_array[j].GetBinCenter(l-1)
			lower_containment_16[j][l] = hist_array[j].Integral(0,l-1)/hist_array[j].Integral()
		if (ratio > 0.05):
			#print hist_array[j].GetBinCenter(l-1)
			lower_bound_05[j][l] = hist_array[j].GetBinCenter(l-1)
			lower_containment_05[j][l] = hist_array[j].Integral(0,l-1)/hist_array[j].Integral()
	#Extract first 2 non-zero values of each array and interpolate to get accurate 0.16 and 0.05 bounds
	if (len(np.trim_zeros(lower_bound_16[j]))>0):
		lower_bound_16_low = np.trim_zeros(lower_bound_16[j])[0]
	else:
		print "No lower bound for:", lamb[j]
	if (len(np.trim_zeros(lower_bound_16[j]))>1):
		lower_bound_16_high= np.trim_zeros(lower_bound_16[j])[1]
	else:
		lower_bound_16_high= 0
	if (len(np.trim_zeros(lower_containment_16[j]))>0):
		lower_containment_16_low = np.trim_zeros(lower_containment_16[j])[0]
	else:
		print "No lower containment for:", lamb[j]
	if (len(np.trim_zeros(lower_containment_16[j]))>1):
		lower_containment_16_high = np.trim_zeros(lower_containment_16[j])[1]
	else:
		lower_containment_16_high = 0
	low_bound_16.append(((lower_bound_16_high-lower_bound_16_low)*0.16/(lower_containment_16_high-lower_containment_16_low))-(lower_containment_16_low*((lower_bound_16_high-lower_bound_16_low)/(lower_containment_16_high-lower_containment_16_low))-lower_bound_16_low))
	if (len(np.trim_zeros(lower_bound_05[j]))>0):
		lower_bound_05_low = np.trim_zeros(lower_bound_05[j])[0]
	else:
		print "No lower bound for:", lamb[j]
	if (len(np.trim_zeros(lower_bound_05[j]))>1):
		lower_bound_05_high= np.trim_zeros(lower_bound_05[j])[1]
	else:
		lower_bound_05_high= 0
	if (len(np.trim_zeros(lower_containment_05[j]))>0):
		lower_containment_05_low = np.trim_zeros(lower_containment_05[j])[0]
	else:
		print "No lower containment for:", lamb[j]
	if (len(np.trim_zeros(lower_containment_05[j]))>1):
		lower_containment_05_high = np.trim_zeros(lower_containment_05[j])[1]
	else:
		lower_containment_05_high = 0
	low_bound_05.append(((lower_bound_05_high-lower_bound_05_low)*0.05/(lower_containment_05_high-lower_containment_05_low))-(lower_containment_05_low*((lower_bound_05_high-lower_bound_05_low)/(lower_containment_05_high-lower_containment_05_low))-lower_bound_05_low))

	#print lower_bound_16_low," ",lower_bound_16_high," ",lower_containment_16_low," ",lower_containment_16_high," ",low_bound_16[j]
	#print lower_bound_05_low," ",lower_bound_05_high," ",lower_containment_05_low," ",lower_containment_05_high," ",low_bound_05[j]
	#Now, higher bounds
	for l in reversed(range(hist_array[j].GetSize())):
		ratio = hist_array[j].Integral(l,hist_array[j].GetSize())/hist_array[j].Integral()
		if (ratio < 0.16):
			higher_bound_16[j][l] = hist_array[j].GetBinCenter(l-1)
			higher_containment_16[j][l] = hist_array[j].Integral(l-1,hist_array[j].GetSize())/hist_array[j].Integral()
		if (ratio < 0.05):
			higher_bound_05[j][l] = hist_array[j].GetBinCenter(l-1)
			higher_containment_05[j][l] = hist_array[j].Integral(l-1,hist_array[j].GetSize())/hist_array[j].Integral()
	if (len(np.trim_zeros(higher_bound_16[j]))>0):
		higher_bound_16_low = np.trim_zeros(higher_bound_16[j])[0]
	else:
		print "No higher bound for:", lamb[j]
	if (len(np.trim_zeros(higher_bound_16[j]))>1):	 
		higher_bound_16_high= np.trim_zeros(higher_bound_16[j])[1]
	else:
		higher_bound_16_high=0
	if (len(np.trim_zeros(higher_containment_16[j]))>0):
		higher_containment_16_low = np.trim_zeros(higher_containment_16[j])[0]
	else:
		print "No higher bound for:", lamb[j]
	if (len(np.trim_zeros(higher_containment_16[j]))>1):
		higher_containment_16_high = np.trim_zeros(higher_containment_16[j])[1]
	else:
		higher_containment_16_high = 0
	high_bound_16.append(((higher_bound_16_high-higher_bound_16_low)*0.16/(higher_containment_16_high-higher_containment_16_low))-(higher_containment_16_low*((higher_bound_16_high-higher_bound_16_low)/(higher_containment_16_high-higher_containment_16_low))-higher_bound_16_low))	

	
	if (len(np.trim_zeros(higher_bound_05[j]))>0):
		higher_bound_05_low = np.trim_zeros(higher_bound_05[j])[0]
	else:
		print "No higher bound for:", lamb[j]
	if (len(np.trim_zeros(higher_bound_05[j]))>1):	 
		higher_bound_05_high= np.trim_zeros(higher_bound_05[j])[1]
	else:
		higher_bound_05_high=0

	if (len(np.trim_zeros(higher_containment_05[j]))>0):
		 higher_containment_05_low = np.trim_zeros(higher_containment_05[j])[0]
	else:
		print "No higher bound for:", lamb[j] 
	if (len(np.trim_zeros(higher_containment_05[j]))>1):
		higher_containment_05_high = np.trim_zeros(higher_containment_05[j])[1]
	else:
		higher_containment_05_high = 0
	high_bound_05.append(((higher_bound_05_high-higher_bound_05_low)*0.05/(higher_containment_05_high-higher_containment_05_low))-(higher_containment_05_low*((higher_bound_05_high-higher_bound_05_low)/(higher_containment_05_high-higher_containment_05_low))-higher_bound_05_low))			

	#print higher_bound_16_low," ",higher_bound_16_high," ",higher_containment_16_low," ",higher_containment_16_high," ",high_bound_16[j],"\n"
	#print higher_bound_05_low," ",higher_bound_05_high," ",higher_containment_05_low," ",higher_containment_05_high," ",high_bound_05[j],"\n"
	
	#Plotting Estetics and bound lines
	hist_array[j].SetFillColorAlpha( 38, 0.5 )
	line_low_array.append(TLine(low_bound_16[j],0,low_bound_16[j],hist_array[j].GetBinContent(hist_array[j].GetXaxis().FindBin(low_bound_16[j]))))
	line_high_array.append(TLine(high_bound_16[j],0,high_bound_16[j],hist_array[j].GetBinContent(hist_array[j].GetXaxis().FindBin(high_bound_16[j]))))
	
	#Define regions to plot with tgraphs
	n_bins_region = 2*(hist_array[j].FindBin(high_bound_16[j])-hist_array[j].FindBin(low_bound_16[j]))
	#print n_bins_region/2
	graph=TGraph(2*n_bins_region)
	graph_array.append(graph)
	for l in range(0,n_bins_region/2):
		graph_array[j].SetPoint(2*l,hist_array[j].GetBinLowEdge(hist_array[j].FindBin(low_bound_16[j])+l),1.*hist_array[j].GetBinContent(hist_array[j].FindBin(low_bound_16[j])+l))

		graph_array[j].SetPoint(2*l+1,hist_array[j].GetBinLowEdge(hist_array[j].FindBin(low_bound_16[j])+l+1),1.*hist_array[j].GetBinContent(hist_array[j].FindBin(low_bound_16[j])+l))	
			
		graph_array[j].SetPoint(2*n_bins_region-2*l-1,hist_array[j].GetBinLowEdge(hist_array[j].FindBin(low_bound_16[j])+l),0.*hist_array[j].GetBinContent(hist_array[j].FindBin(low_bound_16[j])+l))	

		graph_array[j].SetPoint(2*n_bins_region-2*l-2,hist_array[j].GetBinLowEdge(hist_array[j].FindBin(low_bound_16[j])+l+1),0.*hist_array[j].GetBinContent(hist_array[j].FindBin(low_bound_16[j])+l))	


print low_bound_05



##****************PLOTTING**********************##


## Plot bounds in SED space

#lamb_new = np.linspace(1,915,100)
lamb_new= np.logspace(np.log10(1),np.log10(1500),100, endpoint=False)

high_bound_05_spline = spline(lamb,high_bound_05,lamb_new, kind='smoothest', order=2)
low_bound_05_spline = spline(lamb,low_bound_05,lamb_new, kind='smoothest', order=2)

high_bound_16_spline = spline(lamb,high_bound_16,lamb_new, kind='smoothest', order=2)
low_bound_16_spline = spline(lamb,low_bound_16,lamb_new, kind='smoothest', order=2)
fig, ax = plt.subplots(1)
ax.set_yscale('log')
ax.set_xscale('log')
ax.fill_between(lamb_new, high_bound_05_spline, low_bound_05_spline, facecolor='red', alpha=0.5, label='Agreement region 95%',interpolate=True)
ax.fill_between(lamb_new, high_bound_16_spline, low_bound_16_spline, facecolor='blue', alpha=0.5, label='Agreement region 68%',interpolate=True)


	##Add some standard EBL model
I_fra =  EBL.readmodel(model = 'franceschini')
z = 0.031
lmu = np.logspace(-0.1,2.8,50)
Int_fra=I_fra.ebl_array(z,lmu)
Int_fra = Int_fra[0,:]
plt.plot(lmu,Int_fra,label = 'Franceschini'.format(z), lw = 2, color='lime')


ax.legend()
plt.ylabel("$\\nu I_{\\nu}$ [nW$\,$m$^{-2}$sr$^{-1}$]",fontsize=15)
plt.xlabel("$\lambda \, [\mu$m]",fontsize=15)
plt.show()



## Plot histograms 


#legend = TLegend(0.2,0.8,0.3,0.9)
legend = TLegend(0.1,0.8,0.5,0.9)


#print I_1
c0c = TCanvas("c0c","c0c",1400, 850)
#c0c = ROOT.TCanvas("canvas")

c0c.Divide(4,3)
#c0c.cd(1)
#hist_array[0].Draw('h')
legend.AddEntry(graph_array[0],"68% containment region")
for i in range(len(lamb)):
	c0c.cd(i+1)
	hist_array[i].GetXaxis().SetTitle("#nu I")
	hist_array[i].GetYaxis().SetTitle("prob")
	#print low_bound_16[i]
	#c0c.Update()

	line_low_array[i].SetLineColor(2)
	line_low_array[i].SetLineStyle(2)
	line_low_array[i].SetLineWidth(2)
	line_high_array[i].SetLineColor(2)
	line_high_array[i].SetLineStyle(2)
	line_high_array[i].SetLineWidth(2)
	hist_array[i].Draw('h')

	graph_array[i].SetFillColorAlpha(2,0.1)	
	graph_array[i].SetFillStyle(3544)	
	graph_array[i].Draw('f')	 
	line_low_array[i].Draw('same') 
	line_high_array[i].Draw('same')
	legend.Draw('same')


plt.show(block=True)

'''
TCanvas* c0c = new TCanvas("c0c","c0c",1400, 850);
c0c->Divide(4,3);
for i in range(len(lamb)):
      c0c->cd(i+1);
      hinten_weight[i]->SetStats(0);
      hinten_weight[i]->Draw();
      char title[25];
      sprintf(title, "vIv @ #lambda=%0.2f", x[i]); 
      hinten_weight[i]->GetXaxis()->SetTitle(title);
      hinten_weight[i]->GetXaxis()->SetTitleSize(0.05);
      hinten_weight[i]->GetXaxis()->SetLabelSize(0.04);
      hinten_weight[i]->GetYaxis()->SetLabelSize(0.04);
    
#hist_5.Draw()


#	print "model: ", model, " chi: ", chi
'''
