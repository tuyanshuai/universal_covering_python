import subprocess
from meshio.write_mfile import *
import os
def show_mesh(F, V):

    write_mfile("tmp.m", F+1, V)
    subprocess.call(["miniMeshViewer.exe", "tmp.m"])
    os.remove("tmp.m")



