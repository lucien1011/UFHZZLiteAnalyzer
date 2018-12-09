import ROOT,os
from PyUtils.UFTier2 import t2_prefix
from PyUtils.Shell import makedirs

# ____________________________________________________________________________________________________________________________________ ||
inputDir        = t2_prefix+'/store/user/t2/users/archived/dsperka/Run2/HZZ4l/SubmitArea_13TeV/rootfiles_Data80X_4lskim_M17_Feb21/'
inputTreeName   = "passedEvents"
outputDir       = "/raid/raid7/lucien/Higgs/HZZ4l/NTuple/ZPlusX/WrongFC/20181209/SkimTree_WrongFC_Run2016Data_v1/"

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
            12.,
            120.,
            40.,
            0.35,
            0.35,
            outputDir,
            fileName,
            True,
            )
    ana.setDebugMode(False)
    ana.loop(inputDir+fileName,inputTreeName)
