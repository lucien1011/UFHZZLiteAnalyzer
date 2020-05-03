import ROOT,os
from PyUtils.UFTier2 import t2_prefix
from PyUtils.Shell import makedirs

# ____________________________________________________________________________________________________________________________________ ||
inputDir        = t2_prefix+'/store/user/t2/users/klo/Higgs/DarkZ/NTuples/BkgMC_Run2017/'
inputTreeName   = "Ana/passedEvents"
#outputDir       = "/raid/raid7/lucien/Higgs/Zprime-NTuple/20190510/SkimTree_Zprime_Run2017Data_m4l70/"
outputDir       = "/raid/raid7/kshi/Zprime/20200212_Zto4l/M4ToInf/SkimTree_Run2017_MMM_MC/"

fileNames = [
    'ZZTo4L_13TeV_powheg_pythia8_RunIIFall17MiniAOD-94X_mc2017_realistic_v10_ext1-v1.root',
    #'GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8.root',
    #'GluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8.root',
    #'GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8.root',
    #'GluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8.root',
    #'GluGluToContinToZZTo4e_13TeV_MCFM701_pythia8.root',
    #'GluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8.root',
    #'WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8.root',  
    #'WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8.root',
    #'ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8.root',
    #'ttH_HToZZ_4LFilter_M125_13TeV_powheg2_JHUGenV7011_pythia8.root',
    #'GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8.root',  
    #'VBF_HToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8.root',
    ]

# ____________________________________________________________________________________________________________________________________ ||
ROOT.gSystem.Load("include/LiteZto4lM4ToInfTreeProducer_h.so")

makedirs(outputDir)
for fileName in fileNames:
    ana = ROOT.LiteZto4lM4ToInfTreeProducer(
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
            fileName,
            )
    ana.loop(inputDir+fileName,inputTreeName)
