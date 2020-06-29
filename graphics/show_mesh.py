import subprocess
from meshio.write_mfile import *
import os
def show_mesh(F, V, C= np.zeros(0)):

    write_mfile("tmp.m", F+1, V, C)
    subprocess.call(["miniMeshViewer.exe", "tmp.m"])
    os.remove("tmp.m")



