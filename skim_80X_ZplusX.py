import os, time

# ____________________________________________________________________________________________________________ ||
mcInfos = [
        [
            'root://cmsio5.rc.ufl.edu//store/user/t2/users/dsperka/Run2/HZZ4l/SubmitArea_13TeV/rootfiles_MC80X_2lskim_M17_Aug18/',
            [   
                "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv2",
                "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv2",
            ],
        ],

    ]

dataInfos = [
        [   
            'root://cmsio5.rc.ufl.edu//store/user/t2/users/dsperka/Run2/HZZ4l/SubmitArea_13TeV/rootfiles_Data80X_2lskim_M17_Feb21/',
            ['SingleElectron_Run2016-03Feb2017',],
        ],
        [   
            'root://cmsio5.rc.ufl.edu//store/user/t2/users/dsperka/Run2/HZZ4l/SubmitArea_13TeV/rootfiles_Data80X_2lskim_M17_Feb21/',
            ['SingleMuon_Run2016-03Feb2017',],
        ],
        [   
            'root://cmsio5.rc.ufl.edu//store/user/t2/users/dsperka/Run2/HZZ4l/SubmitArea_13TeV/rootfiles_Data80X_2lskim_M17_Feb21/',
            ['DoubleMuon_Run2016-03Feb2017',],
        ],
        [   
            'root://cmsio5.rc.ufl.edu//store/user/t2/users/dsperka/Run2/HZZ4l/SubmitArea_13TeV/rootfiles_Data80X_2lskim_M17_Feb21/',
            ['MuonEG_Run2016-03Feb2017',],
        ],
        [   
            'root://cmsio5.rc.ufl.edu//store/user/t2/users/dsperka/Run2/HZZ4l/SubmitArea_13TeV/rootfiles_Data80X_2lskim_M17_Feb21/',
            ['DoubleEG_Run2016-03Feb2017',],
        ],    
          
    ]

# ____________________________________________________________________________________________________________ ||
# Configuration
infos            = dataInfos
outputDir       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180730/SkimTree_Run2016Data/DarkZ_Z1LSelection/"

if not os.path.exists(os.path.abspath(outputDir)):
    os.makedirs(outputDir)

njobs = 1
if njobs > 6: raise RuntimeError, "Too many resources required"
for inputDir,inputSamples in infos:
    for job in range(1,njobs+1):
        for sample in inputSamples:
            #cmd = './ZX_Ana.exe '+inputDir+'/'+sample+' '+outputDir+sample+' 0 '+str(job)+' '+str(njobs)
            cmd = './ZX_Ana.exe '+inputDir+'/'+sample+' '+outputDir+sample+' 1 '+str(job)+' '+str(njobs)
            print cmd
            os.system(cmd)

