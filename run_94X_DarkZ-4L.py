import os, time

# ____________________________________________________________________________________________________________ ||
# Background sample
#dirMC = '/raid/raid9/ahmad/RUN2_Analyzer/CMSSW_8_0_26_patch1/src/liteUFHZZ4LAnalyzer/Ntuples_Feb2018/'
dirMC = '/raid/raid9/ahmad/RUN2_Analyzer/v2/CMSSW_8_0_26_patch1/src/liteUFHZZ4LAnalyzer/Ntuples_Input/2017/'
#dirMC = '/raid/raid9/ahmad/RUN2_Analyzer/v2/CMSSW_8_0_26_patch1/src/liteUFHZZ4LAnalyzer/Ntuples_Input/2016/'
#dirMC = '/raid/raid9/ahmad/RUN2_Analyzer/v2/CMSSW_8_0_26_patch1/src/liteUFHZZ4LAnalyzer/Ntuples_Input/2016/official_signal_sample/'

samplesMC  = [
    ###2016######
    #"ZZTo4L_13TeV_powheg_pythia8_RunIISummer16MiniAODv2",
    #"ZZTo4L_13TeV_powheg_pythia8_ext1_RunIISummer16MiniAODv2"
    
    ###2017###
    #"GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8",
    #"ZZTo4L_13TeV_powheg_pythia8_RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v2",
    "ZZTo4L_13TeV_powheg_pythia8_RunIIFall17MiniAOD-94X_mc2017_realistic_v10_ext1-v1",
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
inputDir        = dirSig
inputSamples    = samplesSig
#inputDir        = dirMC
#inputSamples    = samplesMC
#outputDir       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180619/SkimTree_v1/"
outputDir       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180619/SkimTree_DarkZ_MZd35_v1/"

njobs = 5
if njobs > 6: raise RuntimeError, "Too many resources required"
for job in range(1,njobs+1):
  for sample in inputSamples:
    #cmd = 'nohup ./ZZ4L_Ana.exe '+inputDir+'/'+sample+' '+outputDir+sample+' 0 '+str(job)+' '+str(njobs)+' >& Dump/'+sample+'_'+str(job)+'.log &'
    #cmd = './ZZ4L_Ana.exe '+inputDir+'/'+sample+' '+outputDir+sample+' 0 '+str(job)+' '+str(njobs)+' >& Dump/'+sample+'_'+str(job)+'.log &'
    cmd = './ZZ4L_Ana.exe '+inputDir+'/'+sample+' '+outputDir+sample+' 0 '+str(job)+' '+str(njobs)
    print cmd
    os.system(cmd)

