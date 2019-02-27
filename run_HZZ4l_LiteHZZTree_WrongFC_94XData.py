import ROOT,os
from PyUtils.UFTier2 import t2_prefix
from PyUtils.Shell import makedirs

# ____________________________________________________________________________________________________________________________________ ||
inputDir        = t2_prefix+'/store/user/t2/users/klo/Higgs/DarkZ/NTuples/ZPlusX_Early2017_v1/'
inputTreeName           = "Ana/passedEvents"
outputDir               = "/raid/raid7/lucien/Higgs/HZZ4l/NTuple/ZPlusX/WrongFC/20181209/SkimTree_WrongFC_Run2017Data_v1/"

fileNames = [
        'SingleElectron_Run2017-17Nov2017.root',
        'SingleMuon_Run2017-17Nov2017.root',
        'DoubleMuon_Run2017-17Nov2017.root',
        'DoubleEG_Run2017-17Nov2017.root',
        'MuonEG_Run2017-17Nov2017.root',
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
            9999999.,
            0.35,
            outputDir,
            fileName,
            True,
            )
    ana.setDebugMode(False)
    ana.loop(inputDir+fileName,inputTreeName)
