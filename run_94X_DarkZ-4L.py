import os, time

# ____________________________________________________________________________________________________________ ||
# Background sample
#dirMC = '/raid/raid9/ahmad/RUN2_Analyzer/CMSSW_8_0_26_patch1/src/liteUFHZZ4LAnalyzer/Ntuples_Feb2018/'
#dirMC = '/raid/raid9/ahmad/RUN2_Analyzer/v2/CMSSW_8_0_26_patch1/src/liteUFHZZ4LAnalyzer/Ntuples_Input/2017/'
#dirMC = '/raid/raid9/ahmad/RUN2_Analyzer/v2/CMSSW_8_0_26_patch1/src/liteUFHZZ4LAnalyzer/Ntuples_Input/2016/'
#dirMC = '/raid/raid9/ahmad/RUN2_Analyzer/v2/CMSSW_8_0_26_patch1/src/liteUFHZZ4LAnalyzer/Ntuples_Input/2016/official_signal_sample/'
dirMC = 'root://cmsio5.rc.ufl.edu//store/user/klo/DarkPhoton_Moriond17_NTuple/BkgMC_Run2017/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/crab_GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/180702_104442/0000/'
#dirMC = 'root://cmsio5.rc.ufl.edu//store/user/klo/DarkPhoton_Moriond17_NTuple/BkgMC_Early2017/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/crab_VBF_HToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/180702_114215/0000/'
#dirMC = 'root://cmsio5.rc.ufl.edu//store/user/klo/DarkPhoton_Moriond17_NTuple/BkgMC_Early2017/WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/crab_WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/180702_114340/0000/'
#dirMC = 'root://cmsio5.rc.ufl.edu//store/user/klo/DarkPhoton_Moriond17_NTuple/BkgMC_Early2017/WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/crab_WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/180702_134408/0000/'
#dirMC = 'root://cmsio5.rc.ufl.edu//store/user/klo/DarkPhoton_Moriond17_NTuple/BkgMC_Early2017/ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8/crab_ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8/180702_134649/0000/' 
samplesMC  = [
    ###2016######
    #"ZZTo4L_13TeV_powheg_pythia8_RunIISummer16MiniAODv2",
    #"ZZTo4L_13TeV_powheg_pythia8_ext1_RunIISummer16MiniAODv2"
    
    ###2017###
    #"GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8",
    #"ZZTo4L_13TeV_powheg_pythia8_RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v2",
    #"ZZTo4L_13TeV_powheg_pythia8_RunIIFall17MiniAOD-94X_mc2017_realistic_v10_ext1-v1",
    "GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8_%s"%i for i in range(1,11) if i != 9
    #"VBF_HToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8_%s"%i for i in range(1,7)
    #"WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8_%s"%i for i in range(1,7) if i != 2
    #"WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8_%s"%i for i in range(2,5)
    #"ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8_%s"%i for i in range(1,20)
]

mcInfos = [
        [
            'root://cmsio5.rc.ufl.edu//store/user/klo/DarkPhoton_Moriond17_NTuple/BkgMC_Run2017/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/crab_GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/180702_104442/0000/',
            ["GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8_%s"%i for i in range(1,11) if i != 9],
        ],
        [
            'root://cmsio5.rc.ufl.edu//store/user/klo/DarkPhoton_Moriond17_NTuple/BkgMC_Early2017/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/crab_VBF_HToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/180702_114215/0000/',
            ["VBF_HToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8_%s"%i for i in range(1,7)],
        ],
        [
            'root://cmsio5.rc.ufl.edu//store/user/klo/DarkPhoton_Moriond17_NTuple/BkgMC_Early2017/WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/crab_WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/180702_114340/0000/',
            ["WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8_%s"%i for i in range(1,7) if i != 2],
        ],
        [
            'root://cmsio5.rc.ufl.edu//store/user/klo/DarkPhoton_Moriond17_NTuple/BkgMC_Early2017/WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/crab_WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/180702_134408/0000/',
            ["WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8_%s"%i for i in range(2,5)],
        ],
        [
            'root://cmsio5.rc.ufl.edu//store/user/klo/DarkPhoton_Moriond17_NTuple/BkgMC_Early2017/ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8/crab_ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8/180702_134649/0000/',
            ["ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8_%s"%i for i in range(1,20)],
        ],
        [
            '/raid/raid9/ahmad/RUN2_Analyzer/v2/CMSSW_8_0_26_patch1/src/liteUFHZZ4LAnalyzer/Ntuples_Input/2017/',
            ["ZZTo4L_13TeV_powheg_pythia8_RunIIFall17MiniAOD-94X_mc2017_realistic_v10_ext1-v1",],
        ],
        [
            'root://cmsio5.rc.ufl.edu//store/user/t2/users/dsperka/Run2/HZZ4l/SubmitArea_13TeV/rootfiles_MC80X_4lskim_M17_Feb21/',
            [   
                "GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8_RunIISummer16MiniAODv2",
                "GluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8_RunIISummer16MiniAODv2",
                "GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8_RunIISummer16MiniAODv2",
                "GluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8_RunIISummer16MiniAODv2",  
                "GluGluToContinToZZTo4e_13TeV_MCFM701_pythia8_RunIISummer16MiniAODv2", 
                "GluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8_RunIISummer16MiniAODv2",
            ],
        ],
        #[
        #    'root://cmsio5.rc.ufl.edu//store/user/t2/users/klo/Higgs/DarkZ/NTuples/BkgMC_Early2017_v1/',
        #    [   
        #        "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8",
        #        "DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8",  
        #        "WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8",
        #    ],
        #],

        ]

dataInfos = [
        [   
            'root://cmsio5.rc.ufl.edu//store/user/t2/users/dsperka/Run2/Zprime/2017/rootfiles_Data_Apr16/',
            ['SingleElectron_Run2017',],
        ],
        [   
            'root://cmsio5.rc.ufl.edu//store/user/t2/users/dsperka/Run2/Zprime/2017/rootfiles_Data_Apr16/',
            ['SingleMuon_Run2017',],
        ],
        [   
            'root://cmsio5.rc.ufl.edu//store/user/t2/users/dsperka/Run2/Zprime/2017/rootfiles_Data_Apr6/',
            ['DoubleMuon_Run2017-17Nov2017-v1',],
        ],
        [   
            'root://cmsio5.rc.ufl.edu//store/user/t2/users/dsperka/Run2/Zprime/2017/rootfiles_Data_Apr6/',
            [       
                'DoubleEG_Run2017B-17Nov2017-v1',  
                'DoubleEG_Run2017E-17Nov2017-v1',   
                'DoubleEG_Run2017D-17Nov2017-v1',  
                'DoubleEG_Run2017C-17Nov2017-v1',  
                'DoubleEG_Run2017F-17Nov2017-v1',   
                ],
        ],
        [   
            'root://cmsio5.rc.ufl.edu//store/user/t2/users/dsperka/Run2/Zprime/2017/rootfiles_Data_Apr6/',
            [       
                'MuonEG_Run2017C-17Nov2017-v1',  
                'MuonEG_Run2017F-17Nov2017-v1',    
                'MuonEG_Run2017B-17Nov2017-v1',  
                'MuonEG_Run2017E-17Nov2017-v1',   
                'MuonEG_Run2017D-17Nov2017-v1',  
                ],
        ],        
          
        ]

# ____________________________________________________________________________________________________________ ||
# Signal sample
#dirSig = 'root://cmsio5.rc.ufl.edu//store/user/klo/DarkPhoton_Moriond17_NTuple/ZD_UpTo0j_MZD20_Eps1e-2/ZD_UpTo0j_MZD20_Eps1e-2/crab_ZD_UpTo0j_MZD20_Eps1e-2_klo/180424_132813/0000/'
#dirSig = 'root://cmsio5.rc.ufl.edu//store/user/klo/DarkPhoton_Moriond17_NTuple/ZD_UpTo0j_MZD15_Eps1e-2/ZD_UpTo0j_MZD15_Eps1e-2/crab_ZD_UpTo0j_MZD15_Eps1e-2_klo/180619_210713/0000/'
#dirSig = 'root://cmsio5.rc.ufl.edu//store/user/klo/DarkPhoton_Moriond17_NTuple/ZD_UpTo0j_MZD25_Eps1e-2/ZD_UpTo0j_MZD25_Eps1e-2/crab_ZD_UpTo0j_MZD25_Eps1e-2_klo/180619_212641/0000/'
#dirSig = 'root://cmsio5.rc.ufl.edu//store/user/klo/DarkPhoton_Moriond17_NTuple/ZD_UpTo0j_MZD30_Eps1e-2/ZD_UpTo0j_MZD30_Eps1e-2/crab_ZD_UpTo0j_MZD30_Eps1e-2_klo/180619_212358/0000/'
dirSig = 'root://cmsio5.rc.ufl.edu//store/user/klo/DarkPhoton_Moriond17_NTuple/ZD_UpTo0j_MZD35_Eps1e-2/ZD_UpTo0j_MZD35_Eps1e-2/crab_ZD_UpTo0j_MZD35_Eps1e-2_klo/180619_213523/0000/'
samplesSig = [
        #"ZD_UpTo0j_MZD20_Eps1e-2_klo_%s"%i for i in range(1,78)
        #"ZD_UpTo0j_MZD15_Eps1e-2_klo_%s"%i for i in range(1,15)
        #"ZD_UpTo0j_MZD25_Eps1e-2_klo_%s"%i for i in range(1,17)
        #"ZD_UpTo0j_MZD30_Eps1e-2_klo_%s"%i for i in range(1,17)
        "ZD_UpTo0j_MZD35_Eps1e-2_klo_%s"%i for i in range(1,18)
        ]

# ____________________________________________________________________________________________________________ ||
# Configuration
infos            = mcInfos
#outputDir       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180729/SkimTree_BkgMC/"
#outputDir       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180730/SkimTree_BkgMC/DarkZ/"
#outputDir       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180730/SkimTree_Run2017Data/DarkZ_m4l105To140_mZ140To120_mZ24To120/"
outputDir       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180730/SkimTree_BkgMC/DarkZ_m4l105To140_mZ140To120_mZ24To120/"
#outputDir       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180729/SkimTree_Run2017Data/"

if not os.path.exists(os.path.abspath(outputDir)):
    os.makedirs(outputDir)

njobs = 1
if njobs > 6: raise RuntimeError, "Too many resources required"
for inputDir,inputSamples in infos:
    for job in range(1,njobs+1):
        for sample in inputSamples:
            cmd = './ZZ4L_Ana.exe '+inputDir+'/'+sample+' '+outputDir+sample+' 0 '+str(job)+' '+str(njobs)
            #cmd = './ZZ4L_Ana.exe '+inputDir+'/'+sample+' '+outputDir+sample+' 1 '+str(job)+' '+str(njobs)
            print cmd
            os.system(cmd)

