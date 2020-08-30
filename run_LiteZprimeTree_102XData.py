import ROOT,os
from PyUtils.UFTier2 import t2_prefix
from PyUtils.Shell import makedirs

# ____________________________________________________________________________________________________________________________________ ||
inputDir        = t2_prefix+'/store/user/t2/users/klo/Higgs/HZZ4l/NTuple/Run2/Data2018_102X_M19_3l_2018jets/'
inputTreeName   = "passedEvents"
outputDir       = "/cmsuf/data/store/user/t2/users/klo/Zprime/EXO-18-008/102X_Data_DarkZNTuple/"

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
            12.,
            999999.,
            0.35,
            outputDir,
            fileName,
            False,
            False,
            )
    ana.loop(inputDir+fileName,inputTreeName)
