import os, time

isDarkPhotonReco = True

if isDarkPhotonReco: 
    dataInfos = [
        [
            #"/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180905/SkimTree_Upsilon_Run2016Data_v1/",
            #["Data_Run2016-03Feb2017_4l_1"],
            #"/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180910/SkimTree_Upsilon_ZX_Run2016Data/",
            #["Data_Run2016_noDuplicates"],
            "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180914/SkimTree_Upsilon_Run2016Data_v1/",
            ["Data_Run2016-03Feb2017_4l_1"],
        ],
        ]
else:
    dataInfos = [
        [
            "root://cmsio5.rc.ufl.edu//store/user/t2/users/archived/dsperka/Run2/HZZ4l/SubmitArea_13TeV/rootfiles_Data80X_4lskim_M17_Feb21/",
            ["Data_Run2016-03Feb2017_4l"],
        ],
        ]


# ____________________________________________________________________________________________________________ ||
# Configuration
infos            = dataInfos
if isDarkPhotonReco: 
    #outputDir       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180905/SkimTree_Upsilon_Run2016Data_v1_DarkPhotonReco/"
    #outputDir       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180910/SkimTree_Upsilon_ZX_Run2016Data_DarkPhotonReco/"
    outputDir       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180914/SkimTree_Upsilon_Run2016Data_v1_DarkPhotonReco/"
else:
    #outputDir       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180905/SkimTree_Upsilon_Run2016Data_v1/"
    outputDir       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180914/SkimTree_Upsilon_Run2016Data_v1/"

if not os.path.exists(os.path.abspath(outputDir)):
    os.makedirs(outputDir)

njobs = 1
if njobs > 6: raise RuntimeError, "Too many resources required"
for inputDir,inputSamples in infos:
    for job in range(1,njobs+1):
        for sample in inputSamples:
            if isDarkPhotonReco:
                cmd = './DarkZ_Ana.exe '+inputDir+'/'+sample+' '+outputDir+sample+' 1 '+str(job)+' '+str(njobs)
            else:
                cmd = './Upsilon_Ana.exe '+inputDir+'/'+sample+' '+outputDir+sample+' 1 '+str(job)+' '+str(njobs)
            os.system(cmd)

