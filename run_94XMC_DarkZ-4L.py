import os, time

t2_dir = 'root://cmsio5.rc.ufl.edu//store/user/t2/users/klo/Higgs/DarkZ/NTuples/BkgMC_Run2017/'

# ____________________________________________________________________________________________________________ ||
mcInfos = [
        [
            t2_dir,
            [
               #'ZZTo4L_13TeV_powheg_pythia8_RunIIFall17MiniAOD-94X_mc2017_realistic_v10_ext1-v1',
               #'WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8',  
               #'WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8',
               #'ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8',
               #'ttH_HToZZ_4LFilter_M125_13TeV_powheg2_JHUGenV7011_pythia8',
               #'GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8',  
               #'VBF_HToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8',
               'GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8',
               'GluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8',
               'GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8',
               'GluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8',
               'GluGluToContinToZZTo4e_13TeV_MCFM701_pythia8',
               'GluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8',
            ],
        ],
    ]

# ____________________________________________________________________________________________________________ ||
# Configuration
infos           = mcInfos

# DarkZ signal region
outputDir       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20181107/SkimTree_DarkPhoton_Run2017Data_m4l70/"

if not os.path.exists(os.path.abspath(outputDir)):
    os.makedirs(outputDir)

njobs = 1
if njobs > 6: raise RuntimeError, "Too many resources required"
for inputDir,inputSamples in infos:
    for job in range(1,njobs+1):
        for sample in inputSamples:
            cmd = './DarkZ_Ana.exe '+inputDir+'/'+sample+' '+outputDir+sample+' 0 '+str(job)+' '+str(njobs)
            os.system(cmd)

