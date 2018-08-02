import os

inputFile       = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180802/SkimTree_Z1LSelection_test/DYJetsToLL_M50_1"
sumWeightFile   = "root://cmsio5.rc.ufl.edu//store/user/t2/users/klo/Higgs/DarkZ/NTuples/ZPlusX_Early2017_v1/DYJetsToLL_M50"
mode            = "DrawZ1LPlot"
outputFile      = "~/public_html/Higgs/DarkZ/ZX/2018-08-02/LeptonPt"

cmd = " ".join(["./ZX_Draw.exe",mode,inputFile,sumWeightFile,outputFile])
print cmd
os.system(cmd)
