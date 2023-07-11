from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH2F
from ROOT import gROOT, gBenchmark, gRandom, gSystem, Double

hist = TH1F('hist', 'my hist', 20, 1, 15)

for i in range(10):
	hist.Fill(i)

c1 = TCanvas( 'c1', 'Histogram Drawing Options', 200, 10, 700, 900 )
#c1.SetFillColor( 1 )
c1.GetFrame().SetFillColor( 2 )
c1.GetFrame().SetBorderSize( 6 )
c1.GetFrame().SetBorderMode( -1 )
hist.Draw()
c1.SaveAs("filename.png")

#c1.Modified()
#c1.Update()
