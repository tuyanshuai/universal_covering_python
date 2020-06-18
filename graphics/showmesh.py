import subprocess
from meshio.write_mfile import *
import os
def showmesh(F, V):

    write_mfile("tmp.m", F,V)
    subprocess.call(["miniMeshViewer.exe", "tmp.m"])
    os.remove("tmp.m")



