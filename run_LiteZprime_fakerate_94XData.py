import ROOT,os
from PyUtils.UFTier2 import t2_prefix
from PyUtils.Shell import makedirs

# ____________________________________________________________________________________________________________________________________ ||
inputDir        = t2_prefix+'/store/user/t2/users/klo/Higgs/DarkZ/NTuples/ZPlusX_Early2017_v1/'
inputTreeName   = "Ana/passedEvents"
#outputDir       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20181116/SkimTree_DarkPhoton_Run2017Data_m4l70/"
outputDir       = "/raid/raid7/kshi/Zprime/20200212_Zto4l/fakerate/SkimTree_Run2017_MMM_Data/"

fileNames = [
        'WZTo3LNu.root',
        'SingleMuon_Run2017-17Nov2017.root',
        'DoubleMuon_Run2017-17Nov2017.root',
        'DoubleEG_Run2017-17Nov2017.root',
        'MuonEG_Run2017-17Nov2017.root',
        'SingleElectron_Run2017-17Nov2017.root',    
        ]

# ____________________________________________________________________________________________________________________________________ ||
ROOT.gSystem.Load("include/LiteZto4lfakerateTreeProducer_h.so")

makedirs(outputDir)
for fileName in fileNames:
    ana = ROOT.LiteZto4lfakerateTreeProducer(
            9999999.,
            70.,
            120.,
            1.,
            120.,
            40.,
            1.,
            999999.,
            0.35,
            outputDir,
            fileName,
            )
    ana.loop(inputDir+fileName,inputTreeName)
