import ROOT,array

xbins_array = array.array('d',[5.,10.,20.,30.,45.,80.,])
fr_el_br = [0.030,0.021,0.0275,0.035,0.062]
fr_el_ec = [0.075,0.0375,0.06,0.095,0.195]
fr_mu_br = [0.100,0.075,0.055,0.070,0.060]
fr_mu_ec = [0.115,0.090,0.060,0.060,0.060]
nbins = 5

h_fr_el_br = ROOT.TH1D("h1D_FRel_EB","",nbins,xbins_array)
h_fr_el_ec = ROOT.TH1D("h1D_FRel_EE","",nbins,xbins_array)
h_fr_mu_br = ROOT.TH1D("h1D_FRmu_EB","",nbins,xbins_array)
h_fr_mu_ec = ROOT.TH1D("h1D_FRmu_EE","",nbins,xbins_array)

for ibin in range(1,nbins+1):
    h_fr_el_br.SetBinContent(ibin,fr_el_br[ibin-1])
    h_fr_el_ec.SetBinContent(ibin,fr_el_ec[ibin-1])
    h_fr_mu_br.SetBinContent(ibin,fr_mu_br[ibin-1])
    h_fr_mu_ec.SetBinContent(ibin,fr_mu_ec[ibin-1])
muFile = ROOT.TFile("/home/lucien/UF-PyNTupleRunner/DarkZ/Data/FakeRate/fakeRates_mu_v2.root","RECREATE")
muFile.cd()
h_fr_mu_br.Write()
h_fr_mu_ec.Write()
muFile.Close()

elFile = ROOT.TFile("/home/lucien/UF-PyNTupleRunner/DarkZ/Data/FakeRate/fakeRates_el_v2.root","RECREATE")
elFile.cd()
h_fr_el_br.Write()
h_fr_el_ec.Write()
elFile.Close()
