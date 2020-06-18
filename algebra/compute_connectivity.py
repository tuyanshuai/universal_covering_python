"""

function [vvif,nvif,pvif] = compute_connectivity(face)            
fi = face(:,1);
fj = face(:,2);
fk = face(:,3);
ff = (1:size(face,1))';
vvif = sparse([fi;fj;fk],[fj;fk;fi],[ff;ff;ff]);
nvif = sparse([ff;ff;ff],[fi;fj;fk],[fj;fk;fi]);
pvif = sparse([ff;ff;ff],[fj;fk;fi],[fi;fj;fk]);

"""
from scipy import sparse
import numpy as np


def compute_connectivity(face):
    fi = face[:, 0]
    fj = face[:, 1]
    fk = face[:, 2]
    ff = np.arange(face.shape[0])


    I = np.array([fi,fj,fk])
    J = np.array([fj,fk,fi])
    V = np.array([ff,ff,ff])
    vvif = sparse.coo_matrix((V.flatten(), (I.flatten(), J.flatten())))
    vvif = vvif.tolil()


    I = np.array([ff,ff,ff])
    J = np.array([fi,fj,fk])
    V = np.array([fj,fk,fi])
    nvif = sparse.coo_matrix((V.flatten(), (I.flatten(), J.flatten())))
    nvif = nvif.tolil()


    I = np.array([ff,ff,ff])
    J = np.array([fj,fk,fi])
    V = np.array([fi,fj,fk])
    pvif = sparse.coo_matrix((V.flatten(), (I.flatten(), J.flatten())))
    pvif = pvif.tolil()


    return vvif,nvif,pvif