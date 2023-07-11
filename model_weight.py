##NOTE: in order to display the histogram, the script must be run with the -i option
import numpy as np
import matplotlib.pyplot as plt 
import os
from ROOT import TH1F, TH1D
from ROOT import TCanvas, TPad, TFile, TPaveLabel, TPaveText, TFile
from ROOT import gROOT

f = TFile("test.png","RECREATE")

hist_1 = TH1F('hist_1', 'Intensity @ #lambda = 1#mum', 50, 1, 20)
hist_2 = TH1F('hist_1', 'Intensity @ #lambda = 84#mum', 50, 1, 20)
hist_3 = TH1F('hist_1', 'Intensity @ #lambda = 500#mum', 50, 1, 20)
hist_4 = TH1F('hist_1', 'Intensity @ #lambda = 833#mum', 50, 1, 20)

f_modelIntensities=open('/home/mateo/Documents/HAWC/EBL/EBLmodel_likelihood/makeSplines/EBLintensities_knotfree_log10_100sample_FULL.dat', 'r')
d_modelIntensities = f_modelIntensities.readlines()
#reader = csv.reader(f, delimiter=' ', skipinitialspace=True)
#first_row = next(reader)
#num_cols = len(first_row)

#lamb=[1,  1.77828,  3.16228,  5.62341,  10,  17.7828,  31.6228,  56.2341, 100,  177.828,  316.228,  562.341]
nmodels=100
#lamb=[5,  6.05764,  7.339,  8.8914,  10.7722,  13.0508,  15.8114,  19.1559,  23.2079,  28.1171,  34.0646,  41.2702]
lamb=[1,  84.25,  167.5,  250.75,  334,  417.25,  500.5,  583.75,  667,  750.25,  833.5,  916.75]

Intensity=np.empty([nmodels,len(lamb)])
i=-1

## Create array with intensities for each model (first index) and wavelenght (second index) 

for x in d_modelIntensities: 
	#Intensity.append([])
        i=i+1;
	for s in range(len(lamb)):
		I=float(x.split('  ')[s])
		Intensity[i][s]=I
#print Intensity

##create map with sigma value in I - lambda space






## Create array
f_results=open('/home/mateo/Documents/HAWC/EBL/EBLmodel_likelihood/mrk421_gp_SSC_EBL_loglike_Fixed/Model_BICSigmaValues.txt', 'r')
d_results = f_results.readlines()
i=-1
#Intensity_w=[]## weight for each model (1st index) and each lamba (2nd index)
Intensity_w=np.empty([nmodels,len(lamb)])


total_points = Intensity.size
Int_array = np.empty([total_points])
lambda_array = np.empty([total_points])
#sigma_array = np.empty([total_points,total_points])
sigma_array = np.empty([total_points])
k=-1
for x in d_results:
	sigma= float(x.split(' ')[1])
	i=i+1	
	for j in range(len(lamb)):
		k=k+1
                Int_array[k]=Intensity[i][j]
                lambda_array[k]=lamb[j]
		sigma_array[k]=100/sigma
		Intensity_w[i][j]=1/sigma
		#sigma=float(x.split(' ')[2])

	#hist_5.Fill(Intensity[i][0])

#sigma_array = np.empty([total_points,total_points])
#print Intensity_w
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
ax.scatter(lambda_array, Int_array, sigma_array,marker='o')
plt.gray()
#cp = ax.contourf(X, Y, sigma_array)
#fig.colorbar(cp) # Add a colorbar to a plot
#ax.set_title('Filled Contours Plot')

plt.ylabel("$\\nu I_{\\nu}$ [nW$\,$m$^{-2}$sr$^{-1}$]",fontsize=15)
plt.xlabel("$\lambda \, [\mu$m]",fontsize=15)


#plt.savefig('AgreementMap.png')


##countour plot
#X, Y = np.meshgrid(lambda_array, Int_array)
#cp = ax.contourf(X, Y, sigma_array)

#print X[12]," ",Y[12]

#fig.colorbar(cp) # Add a colorbar to a plot
#ax.set_title('Filled Contours Plot')

#plt.show()



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
hist_1.GetYaxis().SetTitle("1/#sigma")

hist_2.SetFillColor( 38 )
hist_2.GetXaxis().SetTitle("#nu I")
hist_2.GetYaxis().SetTitle("1/#sigma")

hist_3.SetFillColor( 38 )
hist_3.GetXaxis().SetTitle("#nu I")
hist_3.GetYaxis().SetTitle("1/#sigma")

hist_4.SetFillColor( 38 )
hist_4.GetXaxis().SetTitle("#nu I")
hist_4.GetYaxis().SetTitle("1/#sigma")


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
    '''
#hist_5.Draw()


#	print "model: ", model, " chi: ", chi

