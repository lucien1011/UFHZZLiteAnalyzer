import ROOT,os
from PyUtils.UFTier2 import t2_prefix
from PyUtils.Shell import makedirs

# ____________________________________________________________________________________________________________________________________ ||
inputDir        = t2_prefix+'/store/user/t2/users/rosedj1/Higgs/HZZ4l/NTuple/Run2/MC2018_M19_Mar5_3l_2018Jets_bestCandLegacy/'
inputTreeName   = "Ana/passedEvents"
outputDir       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20190402/SkimTree_DarkPhoton_Run2018Data_m4l70/"

fileNames = [
    'WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8.root',
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
