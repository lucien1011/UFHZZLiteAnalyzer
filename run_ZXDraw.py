import os

# ________________________________________________________________________________________________________________________________________________ ||
# Z1L
#inputFile       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180802/SkimTree_Z1LSelection_test/DYJetsToLL_M50_1"
#sumWeightFile   = "root://cmsio5.rc.ufl.edu//store/user/t2/users/klo/Higgs/DarkZ/NTuples/ZPlusX_Early2017_v1/DYJetsToLL_M50"
#inputFile       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180803/SkimTree_Data80X_Z1LSelection/Data_Run2017_noDuplicates"
#sumWeightFile   = "root://cmsio5.rc.ufl.edu//store/user/t2/users/klo/Higgs/DarkZ/NTuples/ZPlusX_Early2017_v1/DYJetsToLL_M50"

#inputFile       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180806/SkimTree_Data80X_HIG-16-041-ZXCRSelection/Data_Run2016_noDuplicates"
#inputFile       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180806/SkimTree_Data80X_HIG-16-041-ZXCRSelection_v2/Data_Run2016_noDuplicates"
#inputFile       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180820/SkimTree_DarkPhoton_ZX_Run2016Data_v1/Data_Run2016_noDuplicates"
inputFile       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180823/SkimTree_Data80X_HIG-16-041-ZXCRSelectionWithFlag_v3_liteHZZAna/Data_Run2016_noDuplicates_1"
sumWeightFile   = "root://cmsio5.rc.ufl.edu//store/user/t2/users/klo/Higgs/DarkZ/NTuples/ZPlusX_Early2017_v1/DYJetsToLL_M50"
isData          = "1"

#mode            = "DrawZ1LPlot"
mode            = "MakeFRWeight"
#outputFile      = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180820/SkimTree_DarkPhoton_ZX_Run2016Data_v1/Data_Run2016_noDuplicates_FRWeight"
#outputFile      = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180806/SkimTree_Data80X_HIG-16-041-ZXCRSelection_v2/Data_Run2016_noDuplicates_FRWeight_v5"
outputFile      = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180823/SkimTree_Data80X_HIG-16-041-ZXCRSelectionWithFlag_v3_liteHZZAna/Data_Run2016_noDuplicates_1_FRWeight"

# ________________________________________________________________________________________________________________________________________________ ||
dirname = os.path.dirname(os.path.abspath(outputFile+".root"))
if not os.path.exists(dirname):
    os.makedirs(dirname)
cmd = " ".join(["./ZX_Draw.exe",mode,inputFile,sumWeightFile,outputFile,isData])
os.system(cmd)
