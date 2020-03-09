import ROOT,os
from PyUtils.UFTier2 import t2_prefix
from PyUtils.Shell import makedirs
from Emailer.Utils import sendQuickMail,getTimeStamp

# ____________________________________________________________________________________________________________________________________ ||
inputDir        = t2_prefix+'/store/user/t2/users/archived/dsperka/Run2/HZZ4l/SubmitArea_13TeV/rootfiles_Data80X_2lskim_M17_Feb21/'
inputTreeName   = "Ana/passedEvents"
#outputDir       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20190122/SkimTree_DarkPhoton_Run2016Data_m4l70/"
outputDir       = "/raid/raid7/kshi/Zprime/20200212_Zto4l/promptCR/SkimTree_Run2016_MMM_Data/"

fileNames = [
    'SingleMuon_Run2016-03Feb2017.root',
    'DoubleMuon_Run2016-03Feb2017.root',
    ]

# ____________________________________________________________________________________________________________________________________ ||
ROOT.gSystem.Load("include/LiteZto4llowGevTreeProducer_h.so")

makedirs(outputDir)
for fileName in fileNames:
    ana = ROOT.LiteZto4llowGevTreeProducer(
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
