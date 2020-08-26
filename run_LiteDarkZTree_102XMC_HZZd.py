import ROOT,os
from PyUtils.UFTier2 import t2_prefix
from PyUtils.Shell import makedirs

# ____________________________________________________________________________________________________________________________________ ||
inputDir        = t2_prefix+'/store/user/t2/users/klo/Higgs/HToZZd/102X_MCProd_191127/'
inputTreeName   = "Ana/passedEvents"
outputDir       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20191201/SkimTree_DarkPhoton_Run2018Data_m4l70/"

fileNames = [
        "HToZZdTo4L_M125_MZd10_eps1e-2_13TeV_madgraph_pythia8.root",
        "HToZZdTo4L_M125_MZd25_eps1e-2_13TeV_madgraph_pythia8.root",  
        "HToZZdTo4L_M125_MZd3_eps1e-2_13TeV_madgraph_pythia8.root",
        "HToZZdTo4L_M125_MZd15_eps1e-2_13TeV_madgraph_pythia8.root",  
        "HToZZdTo4L_M125_MZd2_eps1e-2_13TeV_madgraph_pythia8.root",   
        "HToZZdTo4L_M125_MZd4_eps1e-2_13TeV_madgraph_pythia8.root",
        "HToZZdTo4L_M125_MZd1_eps1e-2_13TeV_madgraph_pythia8.root",   
        "HToZZdTo4L_M125_MZd30_eps1e-2_13TeV_madgraph_pythia8.root",  
        "HToZZdTo4L_M125_MZd7_eps1e-2_13TeV_madgraph_pythia8.root",
        "HToZZdTo4L_M125_MZd20_eps1e-2_13TeV_madgraph_pythia8.root",  
        "HToZZdTo4L_M125_MZd35_eps1e-2_13TeV_madgraph_pythia8.root",
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
            40.,
            999999.,
            0.35,
            outputDir,
            fileName,
            )
    ana.loop(inputDir+fileName,inputTreeName)
