import ROOT,os
from PyUtils.UFTier2 import t2_prefix
from PyUtils.Shell import makedirs

# ____________________________________________________________________________________________________________________________________ ||
inputDir        = t2_prefix+'/store/user/t2/users/rosedj1/Higgs/DarkZ/NTuples/ppToZZdTo4l/'
inputTreeName   = "Ana/passedEvents"
outputDir       = "/raid/raid7/rosedj1/Higgs/DarkZ-NTuple/20190328/SkimTree_DarkPhoton_Run2017MC_ppToZZdTo4l/"

fileNames = [
    'ppToZZdTo4l_mZd4GeV.root',
    'ppToZZdTo4l_mZd15GeV.root',
    'ppToZZdTo4l_mZd30GeV.root',
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
