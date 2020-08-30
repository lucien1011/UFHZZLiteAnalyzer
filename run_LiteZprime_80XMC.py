import ROOT,os
from PyUtils.UFTier2 import t2_prefix
from PyUtils.Shell import makedirs

# ____________________________________________________________________________________________________________________________________ ||
inputDir        = t2_prefix+'/store/user/t2/users/klo/Zprime/EXO-18-008/HZZNTuple_Run2016/'
inputTreeName   = "Ana/passedEvents"
#outputDir       = "/raid/raid7/lucien/Higgs/Zprime-NTuple/20200415/SkimTree_Zprime_Run2017Data_m4l70/"
outputDir       = "/cmsuf/data/store/user/t2/users/klo/Zprime/EXO-18-008/80X_MCProd_DarkZNTuple/"

fileNames = [
    'ZZTo4L_13TeV_powheg_pythia8_RunIISummer16MiniAODv2.root',
    'GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8_RunIISummer16MiniAODv2.root',
    "ZpTomumu_M10_13TeV_MadGraph5_pythia8-v4_muahmad-RunIISummer16MiniAODv2.root",  
    "ZpTomumu_M50_13TeV_MadGraph5_pythia8-v4_muahmad-RunIISummer16MiniAODv2.root",
    "ZpTomumu_M1_13TeV_MadGraph5_pythia8-v4_muahmad-RunIISummer16MiniAODv2.root",   
    "ZpTomumu_M5_13TeV_MadGraph5_pythia8-v4_muahmad-RunIISummer16MiniAODv2.root",
    "ZpTomumu_M15_13TeV_MadGraph5_pythia8-v4_muahmad-RunIISummer16MiniAODv2.root",  
    "ZpTomumu_M60_13TeV_MadGraph5_pythia8-v4_muahmad-RunIISummer16MiniAODv2.root",
    "ZpTomumu_M20_13TeV_MadGraph5_pythia8-v4_muahmad-RunIISummer16MiniAODv2.root",  
    "ZpTomumu_M70_13TeV_MadGraph5_pythia8-v4_muahmad-RunIISummer16MiniAODv2.root",
    "ZpTomumu_M30_13TeV_MadGraph5_pythia8-v4_muahmad-RunIISummer16MiniAODv2.root",  
    "ZpTomumu_M80_13TeV_MadGraph5_pythia8-v4_muahmad-RunIISummer16MiniAODv2.root",
    "ZpTomumu_M40_13TeV_MadGraph5_pythia8-v4_muahmad-RunIISummer16MiniAODv2.root",
    ]

# ____________________________________________________________________________________________________________________________________ ||
ROOT.gSystem.Load("include/LiteHZZTreeProducer_h.so")

makedirs(outputDir)
for fileName in fileNames:
    ana = ROOT.LiteHZZTreeProducer(
            9999999.,
            70.,
            120.,
            4.,
            120.,
            12.,
            999999.,
            0.35,
            outputDir,
            fileName,
            False,
            False,
            4.0,
            True,
            )
    ana.loop(inputDir+fileName,inputTreeName)
