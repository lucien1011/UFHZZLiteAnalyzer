import ROOT,os
from PyUtils.UFTier2 import t2_prefix
from PyUtils.Shell import makedirs

# ____________________________________________________________________________________________________________________________________ ||
inputDir        = t2_prefix+'/store/user/t2/users/klo/Higgs/DarkZ/NTuples/Data_Run2017/'
inputTreeName   = "Ana/passedEvents"
outputDir       = "/cmsuf/data/store/user/t2/users/klo/Zprime/EXO-18-008/94X_Data_DarkZNTuple/"

fileNames = [
    'SingleMuon_Run2017-17Nov2017-v1.root',
    'DoubleMuon_Run2017-17Nov2017-v1.root',
    #'MuonEG-DoubleEG-SingleElectron_Run2017-17Nov2017-v1.root',
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
            )
    ana.loop(inputDir+fileName,inputTreeName)
