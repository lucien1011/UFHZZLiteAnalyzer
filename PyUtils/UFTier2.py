import subprocess,os

#t2_prefix = "root://cmsio5.rc.ufl.edu:1094/"
t2_prefix = "root://cmsio5.rc.ufl.edu/"

def isRootFile(fileName):
    return fileName.endswith(".root")

def skipTrivial(line):
    return bool(line) and line.split()[-1] != "." and line.split()[-1] != ".." and not "GridFTP Server" in line and 'logged in' not in line

def listdir_uberftp(path,selection=isRootFile):
    if not "ufhpc" in os.environ["HOSTNAME"]:
        cmd = ["uberftp", "cmsio.rc.ufl.edu", "ls %s"%path]
        output = subprocess.Popen(cmd,stdout=subprocess.PIPE).communicate()[0]
        return [l.split()[-1] for l in output.split('\r\n') if skipTrivial(l) and selection(l)]
    else:
        return [l for l in os.listdir(path) if selection(l)]

##____________________________________________________________________________||
class FileInfo(object):
    def __init__(self, path,inUFTier2=False):
        self.path = path
        self.name = os.path.basename(path).replace(".root","")
        self.inUFTier2 = inUFTier2

    def file_path(self):
        if self.inUFTier2:
            return self.uberftp_path()
        else:
            return self.path

    def uberftp_path(self):
        if 'ufhpc' in os.environ['HOSTNAME']:
            return self.path
        else:
            return t2_prefix+self.path.replace("/cmsuf/data","")
