import ROOT,os
from PyUtils.UFTier2 import t2_prefix
from PyUtils.Shell import makedirs
from Emailer.Utils import sendQuickMail,getTimeStamp

# ____________________________________________________________________________________________________________________________________ ||
inputDir        = t2_prefix+'/store/user/t2/users/klo/Higgs/HZZ4l/NTuple/Run2/Data_80XM17_FebCombined/'
inputTreeName   = "passedEvents"
#outputDir       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20190307/SkimTree_DarkPhoton_WrongFC_Run2016Data_m4l70/"
outputDir       = "/raid/raid7/kshi/Zprime/20190718/SkimTree_Run2016_Data/"

fileNames = [
    "Data_Run2016-03Feb2017_4l.root",
    ]

# ____________________________________________________________________________________________________________________________________ ||
ROOT.gSystem.Load("include/LiteHZZTreeProducer_fakerate_h.so")

makedirs(outputDir)
for fileName in fileNames:
    ana = ROOT.LiteHZZTreeProducer_fakerate(
            9999999.,
            70.,
            120.,
            4.,
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
