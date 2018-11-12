import os, time

# ____________________________________________________________________________________________________________ ||
dataInfos = [
        [   
            'root://cmsio5.rc.ufl.edu//store/user/t2/users/klo/Higgs/DarkZ/NTuples/ZPlusX_Early2017_v1/',
            [
                'SingleElectron_Run2017-17Nov2017',
                'SingleMuon_Run2017-17Nov2017',
                'DoubleMuon_Run2017-17Nov2017',
                'DoubleEG_Run2017-17Nov2017',
                'MuonEG_Run2017-17Nov2017',
                ],
        ],
    ]

# ____________________________________________________________________________________________________________ ||
# Configuration
infos = dataInfos
outputDir       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20181107/SkimTree_DarkPhoton_ZX_Run2016Data_m4l70/"

if not os.path.exists(os.path.abspath(outputDir)):
    os.makedirs(outputDir)

njobs = 1
if njobs > 6: raise RuntimeError, "Too many resources required"
for inputDir,inputSamples in infos:
    for job in range(1,njobs+1):
        for sample in inputSamples:
            cmd = './DarkZ_Ana.exe '+inputDir+'/'+sample+' '+outputDir+sample+' 1 '+str(job)+' '+str(njobs)
            os.system(cmd)

