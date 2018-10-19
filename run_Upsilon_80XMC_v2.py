import os, time

isDarkPhotonReco = True

# ____________________________________________________________________________________________________________ ||
#mcInfos = [
#        #[
#        #    #'root://cmsio5.rc.ufl.edu//store/user/t2/users/klo/Higgs/DarkZ/NTuples/ZPlusX_Early2017_v1/',
#        #    'root://cmsio5.rc.ufl.edu//store/user/t2/users/dsperka/Run2/HZZ4l/SubmitArea_13TeV/rootfiles_MC80X_2lskim_M17_Feb21/',
#        #    [ 
#        #        "TTJets_Dilept_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv2",
#        #        "WZTo3LNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv2",
#        #        #"DYJetsToLL_M10To50",
#        #        #"DYJetsToLL_M50",
#        #        #"WZTo3LNu",
#        #        #"TTJets",
#        #    ],
#        #],
#        [
#            #'root://cmsio5.rc.ufl.edu//store/user/t2/users/klo/Higgs/DarkZ/NTuples/ZPlusX_Early2017_v1/',
#            'root://cmsio5.rc.ufl.edu//store/user/t2/users/dsperka/Run2/HZZ4l/SubmitArea_13TeV/rootfiles_MC80X_2lskim_M17_Mar11/',
#            [   
#                "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv2",
#                "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv2",
#            ],
#        ],
#    ]

if not isDarkPhotonReco:
    t2_prefix = 'root://cmsio5.rc.ufl.edu/'
    basedir = t2_prefix+'/store/user/t2/users/archived/dsperka/Run2/HZZ4l/SubmitArea_13TeV/rootfiles_MC80X_4lskim_M17_Feb21/'
    mcInfos = [
#        [
#            basedir,
#            ["GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_RunIISummer16MiniAODv2",],
#        ],
#        [
#            basedir,
#            ["VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_RunIISummer16MiniAODv2",],
#        ],
#        [
#            basedir,
#            ["WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8_RunIISummer16MiniAODv2",],
#        ],
#        [
#            basedir,
#            ["WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8_RunIISummer16MiniAODv2"],
#        ],
#        [
#            basedir,
#            ["ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUgenV6_pythia8_RunIISummer16MiniAODv2"],
#        ],
        [
            basedir,
            #["ZZTo4L_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv2",],
            ["ZZTo4L_13TeV_powheg_pythia8_RunIISummer16MiniAODv2",],
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
else:
    #basedir = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180905/SkimTree_Upsilon_Run2016Data_v1/"
    basedir = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180914/SkimTree_Upsilon_Run2016Data_v1/"
    mcInfos = [
        #[
        #    basedir,
            #["ZZTo4L_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv2",],
            #["ZZTo4L_13TeV_powheg_pythia8_RunIISummer16MiniAODv2_1",],
        #],
        [
            basedir,
            [   
                #"GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8_RunIISummer16MiniAODv2_1",
                #"GluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8_RunIISummer16MiniAODv2_1",
                #"GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8_RunIISummer16MiniAODv2_1",
                "GluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8_RunIISummer16MiniAODv2_1",  
                "GluGluToContinToZZTo4e_13TeV_MCFM701_pythia8_RunIISummer16MiniAODv2_1", 
                "GluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8_RunIISummer16MiniAODv2_1",
            ],
        ],
        ]



# ____________________________________________________________________________________________________________ ||
# Configuration
infos            = mcInfos

if isDarkPhotonReco: 
    #outputDir       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180905/SkimTree_Upsilon_Run2016Data_v1_DarkPhotonReco/"
    outputDir       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180914/SkimTree_Upsilon_Run2016Data_v1_DarkPhotonReco/"
else:
    outputDir       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180914/SkimTree_Upsilon_Run2016Data_v1/"

if not os.path.exists(os.path.abspath(outputDir)):
    os.makedirs(outputDir)

njobs = 1
if njobs > 6: raise RuntimeError, "Too many resources required"
for inputDir,inputSamples in infos:
    for job in range(1,njobs+1):
        for sample in inputSamples:
            if not isDarkPhotonReco:
                cmd = './Upsilon_Ana.exe '+inputDir+'/'+sample+' '+outputDir+sample+' 0 '+str(job)+' '+str(njobs)
            else:
                cmd = './DarkZ_Ana.exe '+inputDir+'/'+sample+' '+outputDir+sample+' 0 '+str(job)+' '+str(njobs)
            os.system(cmd)

