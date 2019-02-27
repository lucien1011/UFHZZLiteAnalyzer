import ROOT,os
from PyUtils.UFTier2 import t2_prefix
from PyUtils.Shell import makedirs

# ____________________________________________________________________________________________________________________________________ ||
inputDir_2lskim         = t2_prefix+'/store/user/t2/users/klo/Higgs/DarkZ/NTuples/ZPlusX_Early2017_v1/'
inputDir_4lskim         = t2_prefix+'/store/user/t2/users/klo/Higgs/DarkZ/NTuples/BkgMC_Run2017/'
inputTreeName           = "Ana/passedEvents"
outputDir               = "/raid/raid7/lucien/Higgs/HZZ4l/NTuple/ZPlusX/WrongFC/20181209/SkimTree_WrongFC_Run2017Data_v1/"

fileNames = [
        inputDir_4lskim+'ZZTo4L_13TeV_powheg_pythia8_RunIIFall17MiniAOD-94X_mc2017_realistic_v10_ext1-v1.root',
        inputDir_2lskim+"TTJets.root",
        inputDir_2lskim+"WZTo3LNu.root",
        inputDir_2lskim+"DYJetsToLL_M10To50.root",
        inputDir_2lskim+"DYJetsToLL_M50.root",
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
            9999999.,
            0.35,
            outputDir,
            os.path.basename(fileName),
            True,
            )
    ana.setDebugMode(False)
    ana.loop(fileName,inputTreeName)
