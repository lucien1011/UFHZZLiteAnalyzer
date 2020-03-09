import ROOT,os
from PyUtils.UFTier2 import t2_prefix
from PyUtils.Shell import makedirs
from Emailer.Utils import sendQuickMail,getTimeStamp

# ____________________________________________________________________________________________________________________________________ ||
inputDir        = t2_prefix+'/store/user/t2/users/klo/Higgs/HZZ4l/NTuple/Run2/Data_80X_2lskim_M17_Feb02/'
inputTreeName   = "Ana/passedEvents"
#outputDir       = "/raid/raid7/kshi/Zprime/20190724/SkimTree_Run2016_signalregion_Data/"
outputDir       = "/raid/raid7/kshi/Zprime/20190827/SkimTree_Run2016_MMM_Data/"

fileNames = [
    "DoubleEG.root",
    "DoubleMuon.root",
    "MuonEG.root",
    "SingleElectron.root",
    "SingleMuon.root",
    ]

# ____________________________________________________________________________________________________________________________________ ||
ROOT.gSystem.Load("include/LiteWto3lMMMTreeProducer_Data_h.so")

makedirs(outputDir)
for fileName in fileNames:
    ana = ROOT.LiteWto3lMMMTreeProducer_Data(
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
