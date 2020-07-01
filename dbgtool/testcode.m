% addpath(genpath('C:/MATLABPackage/'))
% [F,V]=read_off('../data/face.off');
% uv =  disk_conformal_map(F, V);
% % dbg(uv)
% 
% 
% mu = compute_bc(F, uv, V);
% hist(abs(mu),100)


addpath(genpath('C:/MATLABPackage/'))
[F,V]=read_off('../data/maxplanck.nf25k.off');
uvw =  spherical_conformal_map(F, V);
% dbg(uv)
uv(:,1) = uvw(:,1)./(1-uvw(:,3));
uv(:,2) = uvw(:,2)./(1-uvw(:,3));

mu = compute_bc(F, uv, V);
hist(abs(mu),100)