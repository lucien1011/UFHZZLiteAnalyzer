import ROOT,os
from PyUtils.UFTier2 import t2_prefix
from PyUtils.Shell import makedirs

# ____________________________________________________________________________________________________________________________________ ||
inputDir        = t2_prefix+'/store/user/t2/users/klo/Higgs/HZZ4l/NTuple/Run2/Data2018_M19_Feb25_CJSLTEvents_RunABCD_bestCandMela_v2/'
inputTreeName   = "passedEvents"
outputDir       = "/raid/raid9/ahmad/DARK_Photon/CMSSW_8_1_0/src/DarkZ-NTuple/20180304/HZZ_Data_2018/"

fileNames = [
    'Data_Run2018_noDuplicates.root'
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
