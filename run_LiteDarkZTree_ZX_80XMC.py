import ROOT,os
from PyUtils.UFTier2 import t2_prefix
from PyUtils.Shell import makedirs
#from Emailer.Utils import sendQuickMail,getTimeStamp

# ____________________________________________________________________________________________________________________________________ ||
bkgTreeDirT2_Feb21      = t2_prefix+"/store/user/t2/users/klo/Higgs/HZZ4l/NTuple/Run2/MC80X_M17_2l_Feb21/"
bkgTreeDirT2_Aug10      = t2_prefix+"/store/user/t2/users/klo/Higgs/HZZ4l/NTuple/Run2/MC80X_M17_2lskim_Aug10/"
bkgTreeDirT2_June2020   = "/cmsuf/data/store/user/t2/users/klo/Zprime/EXO-18-008/94X_MCProd_200612/"
inputTreeName           = "Ana/passedEvents"
#outputDir               = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20181116/SkimTree_DarkPhoton_ZX_Run2016Data_m4l70/"
outputDir               = "/cmsuf/data/store/user/t2/users/klo/UFNTupleRunner_Storage/kinho.lo/Higgs/DarkZ-NTuple/20181116/SkimTree_DarkPhoton_ZX_Run2016Data_m4l70/"

fileNames = [
        #bkgTreeDirT2_Feb21+"DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root",
        #bkgTreeDirT2_Feb21+"TTJets_Dilept_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8.root",
        #bkgTreeDirT2_Feb21+"WZTo3LNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root",
        #bkgTreeDirT2_Aug10+"DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root",
        #bkgTreeDirT2_June2020+"ZGToLLG_01J_5f_lowMLL_lowGPt_TuneCP5_13TeV-amcatnloFXFX-pythia8_Summer16.root",
        bkgTreeDirT2_June2020+"ZZTo4L_M-1toInf_13TeV_powheg_pythia8_Summer16.root",
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
            0.35,
            0.35,
            outputDir,
            os.path.basename(fileName),
            True,
            )
    ana.setDebugMode(False)
    ana.loop(fileName,inputTreeName)
#sendQuickMail(
#            ["klo@cern.ch",],
#            "UFHZZLiteAnalyzer finished processing ("+getTimeStamp()+") ",
#            "\n".join([
#                "Input directory: "+", ".join([bkgTreeDirT2_Feb21,bkgTreeDirT2_Aug10,]),
#                "Output directory: "+outputDir,
#                "File: "+", ".join(fileNames),
#                ]
#                ),
#            )
