import os

def makedirs(dir):
    if not os.path.exists(dir):
        os.makedirs(os.path.abspath(dir))
