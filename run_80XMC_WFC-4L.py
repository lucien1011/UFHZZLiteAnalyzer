import os, time

# ____________________________________________________________________________________________________________ ||
mcInfos = [
        [
            #'root://cmsio5.rc.ufl.edu//store/user/t2/users/klo/Higgs/DarkZ/NTuples/ZPlusX_Early2017_v1/',
            'root://cmsio5.rc.ufl.edu//store/user/t2/users/archived/dsperka/Run2/HZZ4l/SubmitArea_13TeV/rootfiles_MC80X_2lskim_M17_Feb21/',
            [ 
                "TTJets_Dilept_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv2",
                "WZTo3LNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv2",
                #"DYJetsToLL_M10To50",
                #"DYJetsToLL_M50",
                #"WZTo3LNu",
                #"TTJets",
            ],
        ],
        [
            #'root://cmsio5.rc.ufl.edu//store/user/t2/users/klo/Higgs/DarkZ/NTuples/ZPlusX_Early2017_v1/',
            'root://cmsio5.rc.ufl.edu//store/user/t2/users/archived/dsperka/Run2/HZZ4l/SubmitArea_13TeV/rootfiles_MC80X_2lskim_M17_Mar11/',
            [   
                "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv2",
                "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv2",
            ],
        ],
    ]

# ____________________________________________________________________________________________________________ ||
# Configuration
infos            = mcInfos
# HIG-16-041
#outputDir       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180808/SkimTree_HIG-16-041_Run2016Data_v1/"
# DarkZ signal region
#outputDir       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180808/SkimTree_DarkPhoton_Run2016Data_v1/"
# WFC
outputDir       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180920/SkimTree_WFC_Run2016Data_v1/"

if not os.path.exists(os.path.abspath(outputDir)):
    os.makedirs(outputDir)

njobs = 1
if njobs > 6: raise RuntimeError, "Too many resources required"
for inputDir,inputSamples in infos:
    for job in range(1,njobs+1):
        for sample in inputSamples:
            cmd = './WFC_Ana.exe '+inputDir+'/'+sample+' '+outputDir+sample+' 0 '+str(job)+' '+str(njobs)
            os.system(cmd)

