import ROOT,os
from PyUtils.UFTier2 import t2_prefix
from PyUtils.Shell import makedirs
from Emailer.Utils import sendQuickMail,getTimeStamp

# ____________________________________________________________________________________________________________________________________ ||
inputDir        = '/raid/raid7/lucien/Higgs/HZZ4l/NTuple/ZPlusX/ZXCR/20190128_Run2017_ZXCR-Z1LSkim/'
inputTreeName   = "passedEvents"
outputDir       = "/raid/raid7/lucien/Higgs/HZZ4l/NTuple/ZPlusX/ZXCR/20190313_Run2017_ZXCR-Z1LSkim_LiteHZZTree/"

fileNames = [
    'Data_Run2017-17Nov2017_noDuplicates.root',
    'WZTo3LNu.root',
    'ZZTo4L_ext1.root',
    'TTJets.root',
    'DYJetsToLL_M10To50.root',
    'DYJetsToLL_M50.root',
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
            999999.,
            0.35,
            outputDir,
            fileName,
            False,
            )
    ana.setDebugMode(False)
    ana.loop(inputDir+fileName,inputTreeName)
sendQuickMail(
            ["klo@cern.ch",],
            "UFHZZLiteAnalyzer finished processing ("+getTimeStamp()+") ",
            "\n".join([
                "Input directory: "+", ".join([inputDir,]),
                "Output directory: "+outputDir,
                "File: "+", ".join(fileNames),
                ]
                ),
            )
