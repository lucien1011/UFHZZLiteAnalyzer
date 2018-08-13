import os,argparse

parser = argparse.ArgumentParser()
parser.add_argument('--inputDir',action='store')

option = parser.parse_args()

samples = [
        "GluGluHToZZ",
        "VBF",
        "WminusH",
        "WplusH",
        "ZH",
        "ZZTo4L_13TeV_powheg",
        "GluGluToContinToZZTo2e2mu",
        "GluGluToContinToZZTo2e2tau",
        "GluGluToContinToZZTo2mu2tau",
        "GluGluToContinToZZTo4mu",
        "GluGluToContinToZZTo4e",
        "GluGluToContinToZZTo4tau",
        ]

for sample in samples:
    fileNames = [os.path.join(option.inputDir,n) for n in os.listdir(option.inputDir) if n.endswith(".root") and sample in n]
    if not os.path.exists(os.path.join(option.inputDir,sample)):
        os.makedirs(os.path.join(option.inputDir,sample))
    cmd = 'hadd -f '+os.path.join(option.inputDir,sample,"HaddTree.root")+" "+" ".join(fileNames)
    os.system(cmd)
