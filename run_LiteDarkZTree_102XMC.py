import ROOT,os
from PyUtils.UFTier2 import t2_prefix
from PyUtils.Shell import makedirs

# ____________________________________________________________________________________________________________________________________ ||
inputDir        = t2_prefix+'/store/user/t2/users/klo/Higgs/HZZ4l/NTuple/Run2/MC2018_M19_Feb28_NewSF_2017Jets_bestCandLegacy/'
inputTreeName   = "Ana/passedEvents"
outputDir       = "/raid/raid9/ahmad/DARK_Photon/CMSSW_8_1_0/src/DarkZ-NTuple/20180304/HZZ_MC_2018/"

fileNames = [
#        "GluGluHToZZTo4L_M125_13TeV_powheg2_minloHJ_NNLOPS_JHUgenV702_pythia8.root", 
        "GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8.root",
        "ZZTo4L_TuneCP5_13TeV_powheg_pythia8.root",
        "VBF_HToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8.root",
        "WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8.root",
        "WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8.root",
        "ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8.root",
        "GluGluToContinToZZTo4e_13TeV_MCFM701_pythia8.root",       
        "GluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8.root",   
        "GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8.root",
        "GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8.root",
        "GluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8.root",   
        "GluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8.root",     
        "ttH_HToZZ_4LFilter_M125_13TeV_powheg2_JHUGenV7011_pythia8.root",   
        ]

# ____________________________________________________________________________________________________________________________________ ||
ROOT.gSystem.Load("include/LiteHZZTreeProducer_h.so")

makedirs(outputDir)
for fileName in fileNames:
    ana = ROOT.LiteHZZTreeProducer(outputDir,fileName)
    ana.loop(inputDir+fileName,inputTreeName)
