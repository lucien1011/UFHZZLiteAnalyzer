import ROOT,os
from PyUtils.UFTier2 import t2_prefix
from PyUtils.Shell import makedirs
from Emailer.Utils import sendQuickMail,getTimeStamp

# ____________________________________________________________________________________________________________________________________ ||
inputDir        = t2_prefix+'/store/user/t2/users/klo/Higgs/HZZ4l/NTuple/Run2/Data_80XM17_FebCombined/'
inputTreeName   = "passedEvents"
#outputDir       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20190122/SkimTree_DarkPhoton_Run2016Data_m4l70/"
outputDir       = "/raid/raid7/kshi/Zprime/20200212_Zto4l/SkimTree_Run2016_MMM_Data/"

fileNames = [
    "Data_Run2016-03Feb2017_4l.root",
    ]

# ____________________________________________________________________________________________________________________________________ ||
ROOT.gSystem.Load("include/LiteHZZTreeProducer_h.so")

makedirs(outputDir)
for fileName in fileNames:
    ana = ROOT.LiteHZZTreeProducer(outputDir,fileName)
    ana.loop(inputDir+fileName,inputTreeName)
sendQuickMail(
            ["kshi@cern.ch",],
            "UFHZZLiteAnalyzer finished processing ("+getTimeStamp()+") ",
            "\n".join([
                "Input directory: "+inputDir,
                "Output directory: "+outputDir,
                "File: "+", ".join(fileNames),
                ]
                ),
            )
