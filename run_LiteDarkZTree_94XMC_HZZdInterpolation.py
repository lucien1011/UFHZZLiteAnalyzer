import ROOT,os
from PyUtils.UFTier2 import t2_prefix
from PyUtils.Shell import makedirs

# ____________________________________________________________________________________________________________________________________ ||
inputDir        = t2_prefix+'/store/user/t2/users/klo/Higgs/HToZZd/94X_MCProd_191127/'
inputTreeName   = "Ana/passedEvents"
outputDir       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20191204/SkimTree_DarkPhoton_Run2017Data_m4l70/"

fileNames = [
        "HToZZdTo4L_M125_MZd3_eps1e-2_13TeV_madgraph_pythia8.root",
        "HToZZdTo4L_M125_MZd2_eps1e-2_13TeV_madgraph_pythia8.root",   
        "HToZZdTo4L_M125_MZd4_eps1e-2_13TeV_madgraph_pythia8.root",
        "HToZZdTo4L_M125_MZd1_eps1e-2_13TeV_madgraph_pythia8.root",   
        ]


# ____________________________________________________________________________________________________________________________________ ||
ROOT.gSystem.Load("include/LiteHZZTreeProducer_h.so")

makedirs(outputDir)
for fileName in fileNames:
    ana = ROOT.LiteHZZTreeProducer(
            9999999.,
            70.,
            120.,
            0.1,
            120.,
            40.,
            999999.,
            0.35,
            outputDir,
            fileName,
            False,
            True,
            0.1,
            )
    ana.loop(inputDir+fileName,inputTreeName)
