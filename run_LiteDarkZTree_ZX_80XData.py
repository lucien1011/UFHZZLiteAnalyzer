import ROOT,os
from PyUtils.UFTier2 import t2_prefix
from PyUtils.Shell import makedirs

# ____________________________________________________________________________________________________________________________________ ||
inputDir        = t2_prefix+'/store/user/t2/users/archived/dsperka/Run2/HZZ4l/SubmitArea_13TeV/rootfiles_Data80X_2lskim_M17_Feb21/'
inputTreeName   = "Ana/passedEvents"
outputDir       = "/raid/raid7/kshi/Zprime/20190718/SkimTree_Run2016_Data/"

fileNames = [
    'SingleElectron_Run2016-03Feb2017.root',
    'SingleMuon_Run2016-03Feb2017.root',
    'DoubleMuon_Run2016-03Feb2017.root',
    'MuonEG_Run2016-03Feb2017.root',
    'DoubleEG_Run2016-03Feb2017.root',
    ]

# ____________________________________________________________________________________________________________________________________ ||
ROOT.gSystem.Load("include/LiteHZZTreeProducer_fakerate_h.so")

makedirs(outputDir)
for fileName in fileNames:
    ana = ROOT.LiteHZZTreeProducer_fakerate(outputDir,fileName)
    ana.loop(inputDir+fileName,inputTreeName)
