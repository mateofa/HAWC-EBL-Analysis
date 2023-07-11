##NOTE: in order to display the histogram, the script must be run with the -i option
import numpy as np
import matplotlib.pyplot as plt 
import os
from ROOT import TH1F, TH1D
from ROOT import TCanvas, TPad, TFile, TPaveLabel, TPaveText, TFile
from ROOT import gROOT
from scipy.interpolate import spline
from matplotlib.lines import Line2D

# The PyEBL python class contains methods useful for handling EBL SEDs.
import PyEBL
# The PyBlazar python class contains methods useful for handling blazar spectra
# and their associated optical depth calculations, etc.
import PyBlazar
from ebltable.ebl_from_model import EBL 

f = TFile("test.png","RECREATE")

hist_1 = TH1F('hist_1', 'Intensity @ #lambda = 1#mum', 50, 1, 32)
hist_2 = TH1F('hist_1', 'Intensity @ #lambda = 84#mum', 50, 1, 32)
hist_3 = TH1F('hist_1', 'Intensity @ #lambda = 500#mum', 50, 1, 32)
hist_4 = TH1F('hist_1', 'Intensity @ #lambda = 833#mum', 50, 1, 32)

namewithdir = '/home/mateo/Documents/HAWC/EBL/EBLmodel_likelihood/makeSplines/EBLintensities_knotfree_log10_100sample_fromBestBIC.dat'
f_modelIntensities=open(namewithdir, 'r')

namewithdir_full = '/home/mateo/Documents/HAWC/EBL/EBLmodel_likelihood/makeSplines/EBLintensities_knotfree_log10_100sample_fromBestBIC.dat'
f_modelIntensities_full=open(namewithdir_full, 'r')

d_modelIntensities = f_modelIntensities.readlines()
d_modelIntensities_full = f_modelIntensities_full.readlines()
#reader = csv.reader(f, delimiter=' ', skipinitialspace=True)
#first_row = next(reader)
#num_cols = len(first_row)

#lamb=[1,  1.77828,  3.16228,  5.62341,  10,  17.7828,  31.6228,  56.2341, 100,  177.828,  316.228,  562.341]
nmodels=100
#lamb=[5,  6.05764,  7.339,  8.8914,  10.7722,  13.0508,  15.8114,  19.1559,  23.2079,  28.1171,  34.0646,  41.2702]
lamb=[1,  84.25,  167.5,  250.75,  334,  417.25,  500.5,  583.75,  667,  750.25,  833.5,  916.75]

Intensity=np.empty([nmodels,len(lamb)])
Intensity_full=np.empty([nmodels,len(lamb)])


## Create array with intensities for each model (first index) and wavelenght (second index) 
i=-1
for x in d_modelIntensities: 
	#Intensity.append([])
        i=i+1;
	for s in range(len(lamb)):
		I=float(x.split('  ')[s])
		Intensity[i][s]=I


i=-1
for x in d_modelIntensities_full: 
	#Intensity.append([])
        i=i+1;
	for s in range(len(lamb)):
		I=float(x.split('  ')[s])
		Intensity_full[i][s]=I
##create map with sigma value in I - lambda space






## Create array
f_results=open('/home/mateo/Documents/HAWC/EBL/EBLmodel_likelihood/mrk421_gp_SSC_ECPL_EBL_loglike_Fixed_ZOOM_high/Model_BIC_Categories_Ordered.txt', 'r')
d_results = f_results.readlines()

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
for x in d_results:
	cate= float(x.split('  ')[1])
	i=i+1
	cate_array[i]=cate
	for j in range(len(lamb)):
		k=k+1
                Int_array[k]=Intensity[i][j]
                lambda_array[k]=lamb[j]
		Intensity_w[i][j]=np.exp(-cate)

'''
for x in d_results:
	sigma= float(x.split(' ')[1])
	print sigma
	for i in range(total_points):
		for j in range(total_points):
		#	print i," ",j," ",sigma
			sigma_array[i][j]=sigma
'''

#print sigma_array[50]

###Test gauss scatter

fig, ax = plt.subplots()
#ax.scatter(lambda_array, Int_array, marker='o')

#log10x = np.linspace(0,2.75,1000)
x = np.linspace(1,916.75,1000)

for i in range(nmodels):
	
	eblModel = PyEBL.ebl(namewithdir, i+1)
	#x = 10.**log10x
	y = eblModel.eblSED(x)
		
	if cate_array[i]==3:
		plt.plot(x,y,linewidth=2.,color='chocolate',alpha=0.1, label="Cat 3")
		for j in range(len(lamb)):	
			ax.scatter(lamb[j], Intensity[i][j], color='chocolate',alpha=0.,marker='o')	

	if cate_array[i]==2:
		plt.plot(x,y,linewidth=2.,color='blue',alpha=0.3, label="Cat 2")
		for j in range(len(lamb)):	
			ax.scatter(lamb[j], Intensity[i][j] , color='blue',alpha=0.,marker='o')

	if cate_array[i]==1:
		plt.plot(x,y,linewidth=2.,color='lime',alpha=0.5, label="Cat 1")
		for j in range(len(lamb)):	
			ax.scatter(lamb[j], Intensity[i][j], color='lime',alpha=0.,marker='o')

#	if cate_array[i]==0:
#		print i
#		plt.plot(x,y,linewidth=2.,color='red',alpha=1, label="Prefered")
#		for j in range(len(lamb)):
#			ax.scatter(lamb[j], Intensity[i][j], color='red',alpha=1,marker='o')



###Plot best fit model from full range
f_best=open('/home/mateo/Documents/HAWC/EBL/EBLmodel_likelihood/mrk421_gp_SSC_ECPL_EBL_loglike_Fixed_AllRange/BestModel.txt', 'r')
d_best = f_best.readlines()
for line_best in d_best:
	best_model= int(line_best.split('  ')[3])

print "Current model with min BIC is: ", best_model
eblModel = PyEBL.ebl(namewithdir_full, best_model+1)
y = eblModel.eblSED(x)
plt.plot(x,y,linewidth=3.,color='red',alpha=1, label="Prefered")


### Get traditional models
I_gil =  EBL.readmodel(model = 'franceschini')
lmu = np.logspace(-0.1,2.8,50)
#low_edge=10**(-0.1)
#high_edge=10**2.8
lmu_lin = np.linspace(10**(-0.1),10**(2.8),50)
#print lmu_lin


z=0.031
Int_gil=I_gil.ebl_array(z,lmu)
Int_gil = Int_gil[0,:]

plt.plot(lmu_lin,Int_gil,label = 'Gilmore'.format(z), lw = 2, color='cyan')


###******************************PLOTTING***************************************************###
#print lambda_array[0]," ",Int_array[0]
plt.ylabel("$\\nu I_{\\nu}$ [nW$\,$m$^{-2}$sr$^{-1}$]",fontsize=15)
plt.xlabel("$\lambda \, [\mu$m]",fontsize=15)

##Set legend
custom_lines = [Line2D([0], [0], color='red', lw=3),
                Line2D([0], [0], color='lime', lw=2),
                Line2D([0], [0], color='blue', lw=2),
                Line2D([0], [0], color='chocolate', lw=2),
                Line2D([0], [0], color='cyan', lw=2)]

#fig, ax = plt.subplots()
#lines = ax.plot(data)
plt.legend(custom_lines, ['Best Fit', '2<$\Delta$BIC<6', '6<$\Delta$BIC<10','10<$\Delta$BIC', 'Reference (Franceschini)'])
#plt.legend(loc='best')

#plt.savefig('AgreementMap.png')
#plt.show()

##countour plot
#X, Y = np.meshgrid(lambda_array, Int_array)
#cp = ax.contourf(X, Y, sigma_array)

#print X[12]," ",Y[12]

#fig.colorbar(cp) # Add a colorbar to a plot
#ax.set_title('Filled Contours Plot')


 


#hist_5.Draw()
#f.Write()
##What you want to plot is Intensity_w vs Intensity (for each lambda) Intensity[model][lambda]


#plot intensity weights for each lamba:
I_1=[]
I_2=[]
I_3=[]
I_4=[]

#hist_5.GetXAxis().SetTitle("bla")
#for i in range(len(lamb)):
for i in range(len(d_results)):
	I_1.append(Intensity[i][0])
	I_2.append(Intensity[i][1])
	I_3.append(Intensity[i][6])
	I_4.append(Intensity[i][10])
	hist_1.Fill(Intensity[i][0],Intensity_w[i][0])
	hist_2.Fill(Intensity[i][1],Intensity_w[i][1])
	hist_3.Fill(Intensity[i][6],Intensity_w[i][6])
	hist_4.Fill(Intensity[i][10],Intensity_w[i][10])

hist_1.SetFillColor( 38 )
hist_1.GetXaxis().SetTitle("#nu I")
hist_1.GetYaxis().SetTitle("exp[-cate]")

hist_2.SetFillColor( 38 )
hist_2.GetXaxis().SetTitle("#nu I")
hist_2.GetYaxis().SetTitle("exp[-cate]")

hist_3.SetFillColor( 38 )
hist_3.GetXaxis().SetTitle("#nu I")
hist_3.GetYaxis().SetTitle("exp[-cate]")

hist_4.SetFillColor( 38 )
hist_4.GetXaxis().SetTitle("#nu I")
hist_4.GetYaxis().SetTitle("exp[-cate]")


#print I_1
c0c = TCanvas("c0c","c0c",1400, 850)
c0c.Divide(2,2)
c0c.cd(1)
hist_1.Draw()

c0c.cd(2)
hist_2.Draw()

c0c.cd(3)
hist_3.Draw()

c0c.cd(4)
hist_4.Draw()

c0c.Update()
plt.show()

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
