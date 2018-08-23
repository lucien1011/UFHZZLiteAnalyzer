import os, time

isMC = False
t2_prefix = 'root://cmsio5.rc.ufl.edu/'
basedir = t2_prefix+'/store/user/t2/users/dsperka/Run2/HZZ4l/SubmitArea_13TeV/rootfiles_MC80X_4lskim_M17_Feb21/'

# ____________________________________________________________________________________________________________ ||
mcInfos = [
        [
            basedir,
            ["GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_RunIISummer16MiniAODv2",],
        ],
        [
            basedir,
            ["VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_RunIISummer16MiniAODv2",],
        ],
        [
            basedir,
            ["WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8_RunIISummer16MiniAODv2",],
        ],
        [
            basedir,
            ["WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8_RunIISummer16MiniAODv2"],
        ],
        [
            basedir,
            ["ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUgenV6_pythia8_RunIISummer16MiniAODv2"],
        ],
        [
            basedir,
            ["ZZTo4L_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv2",],
        ],
        [
            basedir,
            [   
                "GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8_RunIISummer16MiniAODv2",
                "GluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8_RunIISummer16MiniAODv2",
                "GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8_RunIISummer16MiniAODv2",
                "GluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8_RunIISummer16MiniAODv2",  
                "GluGluToContinToZZTo4e_13TeV_MCFM701_pythia8_RunIISummer16MiniAODv2", 
                "GluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8_RunIISummer16MiniAODv2",
            ],
        ],
        ]

dataInfos = [
        [
            "root://cmsio5.rc.ufl.edu//store/user/t2/users/dsperka/Run2/HZZ4l/SubmitArea_13TeV/rootfiles_Data80X_4lskim_M17_Feb21/",
            ["Data_Run2016-03Feb2017_4l"],
        ],
        ]

# ____________________________________________________________________________________________________________ ||
# Configuration
infos            = dataInfos if not isMC else mcInfos
outputDir       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180808/SkimTree_DarkPhoton_Run2016Data_v1/"

if not os.path.exists(os.path.abspath(outputDir)):
    os.makedirs(outputDir)

njobs = 1
if njobs > 6: raise RuntimeError, "Too many resources required"
for inputDir,inputSamples in infos:
    for job in range(1,njobs+1):
        for sample in inputSamples:
            if isMC:
                cmd = './DarkZ_Ana.exe '+inputDir+'/'+sample+' '+outputDir+sample+' 0 '+str(job)+' '+str(njobs)
            else:
                cmd = './DarkZ_Ana.exe '+inputDir+'/'+sample+' '+outputDir+sample+' 1 '+str(job)+' '+str(njobs)
            os.system(cmd)

