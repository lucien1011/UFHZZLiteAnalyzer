import ROOT,os
from PyUtils.UFTier2 import t2_prefix
from PyUtils.Shell import makedirs

# ____________________________________________________________________________________________________________________________________ ||
inputDir        = t2_prefix+'/store/user/t2/users/klo/Higgs/HZZ4l/NTuple/Run2/Data_80XM17_FebCombined/'
inputTreeName   = "passedEvents"
#outputDir       = "/raid/raid7/lucien/Higgs/Zprime-NTuple/20200415/SkimTree_Zprime_Run2017Data_m4l70/"
outputDir       = "/cmsuf/data/store/user/t2/users/klo/Zprime/EXO-18-008/80X_Data_DarkZNTuple/"

fileNames = [
    "Data_Run2016-03Feb2017_4l.root",
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
