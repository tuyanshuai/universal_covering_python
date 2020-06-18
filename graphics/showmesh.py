import subprocess
from meshio.write_mfile import *
def showmesh(F, V):

    write_mfile("tmp.m", F,V)

    subprocess.call(["miniMeshViewer.exe", "tmp.m"])



