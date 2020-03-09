import ROOT,os
from PyUtils.UFTier2 import t2_prefix
from PyUtils.Shell import makedirs

# ____________________________________________________________________________________________________________________________________ ||
inputDir        = t2_prefix+'/store/user/t2/users/klo/Higgs/HZZ4l/NTuple/Run2/MC80X_M17_2l_Feb21/'
inputDir2       = t2_prefix+'/store/user/t2/users/klo/Higgs/HZZ4l/NTuple/Run2/MC80X_M17_2lskim_Aug10/'
inputTreeName   = "Ana/passedEvents"
outputDir       = "/raid/raid7/kshi/Zprime/20200212_Zto4l/promptCR/SkimTree_Run2016_MMM_MC/"

fileNames = [
        inputDir+"DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root",
        inputDir+"TTJets_Dilept_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8.root",
        inputDir+"WZTo3LNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root",
        inputDir2+"DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root",
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
            12.,
            1.,
            999999.,
            0.35,
            outputDir,
            os.path.basename(fileName),
            )
    ana.loop(fileName,inputTreeName)
