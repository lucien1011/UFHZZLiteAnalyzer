import ROOT,os
from PyUtils.UFTier2 import t2_prefix
from PyUtils.Shell import makedirs

# ____________________________________________________________________________________________________________________________________ ||
inputDir_2lskimFeb      = t2_prefix+'/store/user/t2/users/archived/dsperka/Run2/HZZ4l/SubmitArea_13TeV/rootfiles_MC80X_2lskim_M17_Feb21/'
inputDir_2lskimMar      = t2_prefix+'/store/user/t2/users/archived/dsperka/Run2/HZZ4l/SubmitArea_13TeV/rootfiles_MC80X_2lskim_M17_Mar11/'
inputDir_4lskimFeb      = t2_prefix+'/store/user/t2/users/archived/dsperka/Run2/HZZ4l/SubmitArea_13TeV/rootfiles_MC80X_4lskim_M17_Feb21/'
inputTreeName           = "Ana/passedEvents"
outputDir               = "/raid/raid7/lucien/Higgs/HZZ4l/NTuple/ZPlusX/WrongFC/20181209/SkimTree_WrongFC_Run2016Data_v1/"

fileNames = [
        inputDir_2lskimFeb+"TTJets_Dilept_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv2.root",
        inputDir_2lskimFeb+"WZTo3LNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv2.root",
        inputDir_2lskimMar+"DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv2.root",
        inputDir_2lskimMar+"DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv2.root",
        inputDir_4lskimFeb+"ZZTo4L_13TeV_powheg_pythia8_RunIISummer16MiniAODv2.root",
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
            0.35,
            0.35,
            outputDir,
            os.path.basename(fileName),
            True,
            )
    ana.setDebugMode(False)
    ana.loop(fileName,inputTreeName)
