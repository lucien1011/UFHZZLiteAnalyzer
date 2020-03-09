import ROOT,os
from PyUtils.UFTier2 import t2_prefix
from PyUtils.Shell import makedirs

# ____________________________________________________________________________________________________________________________________ ||
inputDir        = t2_prefix+'/store/user/t2/users/klo/Higgs/HZZ4l/NTuple/Run2/MC80X_M17_4l_Feb21/'
inputTreeName   = "Ana/passedEvents"
#outputDir       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20190122/SkimTree_DarkPhoton_Run2016Data_m4l70/"
outputDir       = "/raid/raid7/kshi/Zprime/20200212_Zto4l/SkimTree_Run2016_MMM_MC/"

fileNames = [
        #"GluGluHToZZTo4L_M125_13TeV_powheg2_minloHJ_NNLOPS_JHUgenV702_pythia8.root", 
        #"GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8.root",
        "ZZTo4L_13TeV_powheg_pythia8.root",
        #"VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8.root",
        #"WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8.root",
        #"WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8.root",
        #"ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUgenV6_pythia8.root",
        "GluGluToContinToZZTo4e_13TeV_MCFM701_pythia8.root",       
        "GluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8.root",   
        "GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8.root",
        "GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8.root",
        "GluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8.root",   
        "GluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8.root",     
        #"ttH_HToZZ_4LFilter_M125_13TeV_powheg2_JHUgenV6_pythia8.root",   
        ]

# ____________________________________________________________________________________________________________________________________ ||
ROOT.gSystem.Load("include/LiteHZZTreeProducer_h.so")

makedirs(outputDir)
for fileName in fileNames:
    ana = ROOT.LiteHZZTreeProducer(outputDir,fileName)
    ana.loop(inputDir+fileName,inputTreeName)
