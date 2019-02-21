import ROOT,os,glob
from PyUtils.UFTier2 import t2_prefix,listdir_uberftp
from PyUtils.Shell import makedirs
from Emailer.Core import Emailer

# ____________________________________________________________________________________________________________________________________ ||
#baseInputDir    = '/store/user/t2/users/klo/Higgs/HZZ4l/NTuple/Run2/ZXData_Run2018/'
baseInputDir    = '/raid/raid7/lucien/Higgs/HZZ4l/NTuple/ZPlusX/ZXCR/SkimTree_ZX_Run2018Data_190220/'
inputDir        = t2_prefix+baseInputDir
#inputTreeName   = "Ana/passedEvents"
inputTreeName   = "passedEvents"
outputDir       = "/raid/raid7/lucien/Higgs/HZZ4l/NTuple/ZPlusX/ZXCR/SkimTree_WrongFC_Run2018Data_190221/"
sendMail        = False

#fileNames = listdir_uberftp('/cms/data'+baseInputDir)
fileNames = glob.glob(baseInputDir+"*HZZNTuple*_noDuplicates.root")
fileNames = [os.path.basename(f) for f in fileNames]

# ____________________________________________________________________________________________________________________________________ ||
ROOT.gSystem.Load("include/LiteHZZTreeProducer_h.so")

makedirs(outputDir)
for fileName in fileNames:
    print "="*40
    print "Processing "+fileName
    ana = ROOT.LiteHZZTreeProducer(
            9999999.,
            70.,
            120.,
            12.,
            120.,
            40.,
            9999999.,
            0.35,
            outputDir,
            fileName,
            True,
            )
    ana.setDebugMode(False)
    #ana.loop(inputDir+fileName,inputTreeName)
    ana.loop(baseInputDir+fileName,inputTreeName)
if sendMail:
    emailer = Emailer()
    emailer.sendMail(
            "lucien@newberry.ihepa.ufl.edu",
            ["klo@cern.ch","vukasin.milosevic@cern.ch"],
            "UFHZZLiteAnalyzer finished processing ("+emailer.getTimeStr()+") ",
            "\n".join([
                "Input directory: "+baseInputDir,
                "Output directory: "+outputDir,
                "File: "+",".join(fileNames),
                ]
                ),
            )
