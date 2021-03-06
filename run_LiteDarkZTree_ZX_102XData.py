import ROOT,os
from PyUtils.UFTier2 import t2_prefix
from PyUtils.Shell import makedirs

# ____________________________________________________________________________________________________________________________________ ||
inputDir        = t2_prefix+'/store/user/t2/users/klo/Higgs/HZZ4l/NTuple/Run2/Data2018_102X_M19_3l_2018jets/'
inputTreeName   = "Ana/passedEvents"
outputDir       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20190402/SkimTree_DarkPhoton_Run2018Data_m4l70/"

fileNames = [
    'Data_Run2018A-17Sep2018_noDuplicates.root',
    'Data_Run2018B-17Sep2018_noDuplicates.root',
    'Data_Run2018C-17Sep2018_noDuplicates.root',
    'Data_Run2018D-PromptReco-v2_noDuplicates.root',
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
