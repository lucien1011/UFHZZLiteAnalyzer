import os, time

t2_prefix = 'root://cmsio5.rc.ufl.edu/'
basedir = t2_prefix+'/store/user/t2/users/klo/Higgs/DarkZ/NTuples/SigMC_Run2016_v1/'

# ____________________________________________________________________________________________________________ ||
# Signal sample
mcInfos = [
        [
            basedir,
            [
                "ZD_UpTo0j_MZD4_Eps1e-2_klo",
                "ZD_UpTo0j_MZD7_Eps1e-2_klo",
                "ZD_UpTo0j_MZD10_Eps1e-2_klo",
                "ZD_UpTo0j_MZD15_Eps1e-2_klo",
                "ZD_UpTo0j_MZD20_Eps1e-2_klo",
                "ZD_UpTo0j_MZD25_Eps1e-2_klo",
                "ZD_UpTo0j_MZD30_Eps1e-2_klo",
            ],
        ],
    ]

# ____________________________________________________________________________________________________________ ||
# Configuration
infos            = mcInfos
# HIG-16-041
#outputDir       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180808/SkimTree_HIG-16-041_Run2016Data_v1/"
# DarkZ signal region
outputDir       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20181019/SkimTree_DarkPhoton_Run2017Sig_m4l70/"

if not os.path.exists(os.path.abspath(outputDir)):
    os.makedirs(outputDir)

njobs = 1
if njobs > 6: raise RuntimeError, "Too many resources required"
for inputDir,inputSamples in infos:
    for job in range(1,njobs+1):
        for sample in inputSamples:
            cmd = './DarkZ_Ana.exe '+inputDir+'/'+sample+' '+outputDir+sample+' 0 '+str(job)+' '+str(njobs)
            os.system(cmd)

