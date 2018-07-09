import os, time

# ____________________________________________________________________________________________________________ ||
# Background sample
#dirMC = '/raid/raid9/ahmad/RUN2_Analyzer/CMSSW_8_0_26_patch1/src/liteUFHZZ4LAnalyzer/Ntuples_Feb2018/'
#dirMC = '/raid/raid9/ahmad/RUN2_Analyzer/v2/CMSSW_8_0_26_patch1/src/liteUFHZZ4LAnalyzer/Ntuples_Input/2017/'
#dirMC = '/raid/raid9/ahmad/RUN2_Analyzer/v2/CMSSW_8_0_26_patch1/src/liteUFHZZ4LAnalyzer/Ntuples_Input/2016/'
#dirMC = '/raid/raid9/ahmad/RUN2_Analyzer/v2/CMSSW_8_0_26_patch1/src/liteUFHZZ4LAnalyzer/Ntuples_Input/2016/official_signal_sample/'
#dirMC = 'root://cmsio5.rc.ufl.edu//store/user/klo/DarkPhoton_Moriond17_NTuple/BkgMC_Run2017/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/crab_GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/180702_104442/0000/'
#dirMC = 'root://cmsio5.rc.ufl.edu//store/user/klo/DarkPhoton_Moriond17_NTuple/BkgMC_Early2017/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/crab_VBF_HToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/180702_114215/0000/'
#dirMC = 'root://cmsio5.rc.ufl.edu//store/user/klo/DarkPhoton_Moriond17_NTuple/BkgMC_Early2017/WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/crab_WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/180702_114340/0000/'
#dirMC = 'root://cmsio5.rc.ufl.edu//store/user/klo/DarkPhoton_Moriond17_NTuple/BkgMC_Early2017/WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/crab_WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/180702_134408/0000/'
dirMC = 'root://cmsio5.rc.ufl.edu//store/user/klo/DarkPhoton_Moriond17_NTuple/BkgMC_Early2017/ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8/crab_ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8/180702_134649/0000/' 
samplesMC  = [
    ###2016######
    #"ZZTo4L_13TeV_powheg_pythia8_RunIISummer16MiniAODv2",
    #"ZZTo4L_13TeV_powheg_pythia8_ext1_RunIISummer16MiniAODv2"
    
    ###2017###
    #"GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8",
    #"ZZTo4L_13TeV_powheg_pythia8_RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v2",
    #"ZZTo4L_13TeV_powheg_pythia8_RunIIFall17MiniAOD-94X_mc2017_realistic_v10_ext1-v1",
    #"GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8_%s"%i for i in range(1,11) if i != 9
    #"VBF_HToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8_%s"%i for i in [1,5,6]
    #"WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8_%s"%i for i in range(1,7) if i != 2
    #"WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8_%s"%i for i in range(2,5)
    "ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8_%s"%i for i in range(1,20)
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
# Data sample
dirData = "/raid/raid9/ahmad/RUN2_Analyzer/v2/CMSSW_8_0_26_patch1/src/liteUFHZZ4LAnalyzer/Ntuples_Input/2017/"
samplesData = [
'SingleDoubleMuon_Run2017-17Nov2017-v1_NoDuplicates'
#'Data_Run2016-03Feb2017_4l'
]

# ____________________________________________________________________________________________________________ ||
# Configuration
#inputDir        = dirSig
#inputSamples    = samplesSig
inputDir        = dirMC
inputSamples    = samplesMC
#outputDir       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180619/SkimTree_v1/"
outputDir       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180702/SkimTree_BkgMC/"

njobs = 5
if njobs > 6: raise RuntimeError, "Too many resources required"
for job in range(1,njobs+1):
  for sample in inputSamples:
    #cmd = 'nohup ./ZZ4L_Ana.exe '+inputDir+'/'+sample+' '+outputDir+sample+' 0 '+str(job)+' '+str(njobs)+' >& Dump/'+sample+'_'+str(job)+'.log &'
    #cmd = './ZZ4L_Ana.exe '+inputDir+'/'+sample+' '+outputDir+sample+' 0 '+str(job)+' '+str(njobs)+' >& Dump/'+sample+'_'+str(job)+'.log &'
    cmd = './ZZ4L_Ana.exe '+inputDir+'/'+sample+' '+outputDir+sample+' 0 '+str(job)+' '+str(njobs)
    print cmd
    os.system(cmd)

