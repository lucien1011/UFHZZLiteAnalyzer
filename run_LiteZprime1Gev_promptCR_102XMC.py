import ROOT,os
from PyUtils.UFTier2 import t2_prefix
from PyUtils.Shell import makedirs

# ____________________________________________________________________________________________________________________________________ ||
inputDir        = t2_prefix+'/store/user/t2/users/rosedj1/Higgs/HZZ4l/NTuple/Run2/MC2018_M19_Mar12_4l_2018Jets_JER_bestCandLegacy/'
inputTreeName   = "Ana/passedEvents"
#outputDir       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20190402/SkimTree_DarkPhoton_Run2018Data_m4l70/"
outputDir       = "/raid/raid7/kshi/Zprime/20200212_Zto4l/mllLowGev/SkimTree_Run2018_MMM_MC/"

fileNames = [
                'ZZTo4L_TuneCP5_13TeV_powheg_pythia8.root',
                #'GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8.root',  
                #'VBF_HToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8.root',
                #'WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8.root',
                #'WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8.root',
                #'ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8.root',
                #'ttH_HToZZ_4LFilter_M125_13TeV_powheg2_JHUGenV7011_pythia8.root',  
                'GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8.root',         
                'GluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8.root',         
                'GluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8.root',         
                'GluGluToContinToZZTo4e_13TeV_MCFM701_pythia8.root',             
                'GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8.root',             
                'GluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8.root',         
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
