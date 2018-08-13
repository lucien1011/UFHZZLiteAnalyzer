import os, time

isMC = False

# ____________________________________________________________________________________________________________ ||
mcInfos = [
        [
            #'root://cmsio5.rc.ufl.edu//store/user/t2/users/klo/Higgs/DarkZ/NTuples/ZPlusX_Early2017_v1/',
            'root://cmsio5.rc.ufl.edu//store/user/t2/users/dsperka/Run2/HZZ4l/SubmitArea_13TeV/rootfiles_MC80X_2lskim_M17_Feb21/',
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
            'root://cmsio5.rc.ufl.edu//store/user/t2/users/dsperka/Run2/HZZ4l/SubmitArea_13TeV/rootfiles_MC80X_2lskim_M17_Mar11/',
            [   
                "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv2",
                "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv2",
            ],
        ],
    ]

#dataInfos = [
#        [   
#            'root://cmsio5.rc.ufl.edu//store/user/t2/users/dsperka/Run2/HZZ4l/SubmitArea_13TeV/rootfiles_Data80X_2lskim_M17_Feb21/',
#            ['SingleElectron_Run2016-03Feb2017',],
#        ],
#        [   
#            'root://cmsio5.rc.ufl.edu//store/user/t2/users/dsperka/Run2/HZZ4l/SubmitArea_13TeV/rootfiles_Data80X_2lskim_M17_Feb21/',
#            ['SingleMuon_Run2016-03Feb2017',],
#        ],
#        [   
#            'root://cmsio5.rc.ufl.edu//store/user/t2/users/dsperka/Run2/HZZ4l/SubmitArea_13TeV/rootfiles_Data80X_2lskim_M17_Feb21/',
#            ['DoubleMuon_Run2016-03Feb2017',],
#        ],
#        [   
#            'root://cmsio5.rc.ufl.edu//store/user/t2/users/dsperka/Run2/HZZ4l/SubmitArea_13TeV/rootfiles_Data80X_2lskim_M17_Feb21/',
#            ['MuonEG_Run2016-03Feb2017',],
#        ],
#        [   
#            'root://cmsio5.rc.ufl.edu//store/user/t2/users/dsperka/Run2/HZZ4l/SubmitArea_13TeV/rootfiles_Data80X_2lskim_M17_Feb21/',
#            ['DoubleEG_Run2016-03Feb2017',],
#        ],    
#          
#    ]
dataInfos = [
        [
            "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180808/SkimTree_Data80X_HIG-16-041-ZXCRSelectionWithFlag_v3/",
            ["Data_Run2016_noDuplicates",],
        ],
        ]


# ____________________________________________________________________________________________________________ ||
# Configuration
if not isMC:
    infos = dataInfos
else:
    infos = mcInfos
#outputDir       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180803/SkimTree_Data80X_Z1LSelection/"
#outputDir       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180806/SkimTree_MC80X_ZXCRSelection/"
#outputDir       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180806/SkimTree_Data80X_ZXCRSelection/"
outputDir       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180808/SkimTree_Data80X_HIG-16-041-ZXCRSelectionWithFlag_v3_recalculate/"

if not os.path.exists(os.path.abspath(outputDir)):
    os.makedirs(outputDir)

njobs = 1
if njobs > 6: raise RuntimeError, "Too many resources required"
for inputDir,inputSamples in infos:
    for job in range(1,njobs+1):
        for sample in inputSamples:
            #cmd = './ZX_Ana.exe '+inputDir+'/'+sample+' '+outputDir+sample+' 0 '+str(job)+' '+str(njobs)
            #cmd = './ZX_Ana.exe '+inputDir+'/'+sample+' '+outputDir+sample+' 1 '+str(job)+' '+str(njobs)
            cmd = './DarkZ_Ana.exe '+inputDir+'/'+sample+' '+outputDir+sample+' 1 '+str(job)+' '+str(njobs)
            #cmd = './DarkZ_Ana.exe '+inputDir+'/'+sample+' '+outputDir+sample+' 0 '+str(job)+' '+str(njobs)
            print cmd
            os.system(cmd)

