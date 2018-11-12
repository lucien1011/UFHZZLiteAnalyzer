import os

# ________________________________________________________________________________________________________________________________________________ ||
inputFile       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20181107/SkimTree_DarkPhoton_ZX_Run2017Data_m4l70/Data_Run2017_noDuplicates"
sumWeightFile   = "root://cmsio5.rc.ufl.edu//store/user/t2/users/klo/Higgs/DarkZ/NTuples/ZPlusX_Early2017_v1/DYJetsToLL_M50"
isData          = "1"
inputType       = "liteHZZ"

mode            = "MakeFRWeight"
outputFile      = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20181107/SkimTree_DarkPhoton_ZX_Run2017Data_m4l70/Data_Run2017_noDuplicates_FRWeight"


elFilePath      = "/home/lucien/Higgs/DarkZ/CMSSW_9_4_2/src/liteUFHZZ4LAnalyzer/Data/fakeRate2017.root"
muFilePath      = "/home/lucien/Higgs/DarkZ/CMSSW_9_4_2/src/liteUFHZZ4LAnalyzer/Data/fakeRate2017.root"

# ________________________________________________________________________________________________________________________________________________ ||
dirname = os.path.dirname(os.path.abspath(outputFile+".root"))
if not os.path.exists(dirname):
    os.makedirs(dirname)
cmd = " ".join(["./ZX_Weight.exe",mode,inputFile,sumWeightFile,outputFile,isData,inputType,elFilePath,muFilePath,])
os.system(cmd)
