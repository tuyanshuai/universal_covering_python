addpath(genpath('C:/MATLABPackage/'))
[F,V]=read_off('../data/maxplanck.nf25k.off');
uvw =  spherical_conformal_map(F, V);

plot_mesh(F, uvw)