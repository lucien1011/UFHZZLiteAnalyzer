import os

# ________________________________________________________________________________________________________________________________________________ ||
# Z1L
#inputFile       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180802/SkimTree_Z1LSelection_test/DYJetsToLL_M50_1"
#sumWeightFile   = "root://cmsio5.rc.ufl.edu//store/user/t2/users/klo/Higgs/DarkZ/NTuples/ZPlusX_Early2017_v1/DYJetsToLL_M50"
#inputFile       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180803/SkimTree_Data80X_Z1LSelection/Data_Run2017_noDuplicates"
#sumWeightFile   = "root://cmsio5.rc.ufl.edu//store/user/t2/users/klo/Higgs/DarkZ/NTuples/ZPlusX_Early2017_v1/DYJetsToLL_M50"

#inputFile       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180806/SkimTree_Data80X_HIG-16-041-ZXCRSelection/Data_Run2016_noDuplicates"
#inputFile       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180806/SkimTree_Data80X_HIG-16-041-ZXCRSelection_v2/Data_Run2016_noDuplicates"
inputFile       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180820/SkimTree_DarkPhoton_ZX_Run2016Data_v1/Data_Run2016_noDuplicates"
#inputFile       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180823/SkimTree_Data80X_HIG-16-041-ZXCRSelectionWithFlag_v3_liteHZZAna/Data_Run2016_noDuplicates_1"
#inputFile       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180910/SkimTree_Upsilon_ZX_Run2016Data_DarkPhotonReco/Data_Run2016_noDuplicates_1"
#inputFile       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180914/SkimTree_Upsilon_ZX_Run2016Data_DarkPhotonReco/Data_Run2016_noDuplicates_1"
#inputFile       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180920/SkimTree_WFC_Run2016Data_v1/Data_Run2016-03Feb2017_4l_1" 
#inputFile       = "/raid/raid5/predragm/Run2/HZZ4l/SubmitArea_13TeV/rootfiles_MC80X_2lskim_M17_Feb21/Data_ZX_Run2017-03Feb2017_slimmedZX"
inputFile       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180905/SkimTree_DarkPhoton_ZX_Run2016Data_m4l70/Data_Run2016_noDuplicates"
sumWeightFile   = "root://cmsio5.rc.ufl.edu//store/user/t2/users/klo/Higgs/DarkZ/NTuples/ZPlusX_Early2017_v1/DYJetsToLL_M50"
isData          = "1"
#inputType       = "PedjaSkim"
inputType       = "liteHZZ"

#mode            = "DrawZ1LPlot"
mode            = "MakeFRWeight"
#outputFile      = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180820/SkimTree_DarkPhoton_ZX_Run2016Data_v1/Data_Run2016_noDuplicates_FRWeightCorr"
#outputFile      = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180806/SkimTree_Data80X_HIG-16-041-ZXCRSelection_v2/Data_Run2016_noDuplicates_FRWeight_v5"
#outputFile      = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180823/SkimTree_Data80X_HIG-16-041-ZXCRSelectionWithFlag_v3_liteHZZAna/Data_Run2016_noDuplicates_1_FRWeight"
#outputFile       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180910/SkimTree_Upsilon_ZX_Run2016Data_DarkPhotonReco/Data_Run2016_noDuplicates_1_FRWeight"
#outputFile       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180920/SkimTree_WFC_Run2016Data_v1/Data_Run2016-03Feb2017_4l_1_FRWeight"
#outputFile       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180924/PedjaInput_rootfiles_MC80X_2lskim_M17_Feb21/Data_ZX_Run2017-03Feb2017_slimmedZX_FRWeight"
#outputFile      = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180924/SkimTree_DarkPhoton_ZX_Run2016Data_v1/Data_Run2016_noDuplicates_FRWeightCorr"
outputFile      = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180924/SkimTree_DarkPhoton_ZX_Run2016Data_m4l70/Data_Run2016_noDuplicates_FRWeightCorrMyFR"

# ________________________________________________________________________________________________________________________________________________ ||
dirname = os.path.dirname(os.path.abspath(outputFile+".root"))
if not os.path.exists(dirname):
    os.makedirs(dirname)
cmd = " ".join(["./ZX_Weight.exe",mode,inputFile,sumWeightFile,outputFile,isData,inputType])
os.system(cmd)
