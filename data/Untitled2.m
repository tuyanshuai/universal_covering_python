clear all;
dbstop if error
[F,V]=read_obj('bunnyh.obj');

bd = compute_bd(F);
uv = disk_harmonic_map(F,V);
disk = uv(bd, :);
pd = power_diagram(F,uv);


in2 = isinpolygon(disk,pd.dpe);
% 
% nc = size(uv,1);
% 
nc = size(uv,1);
area = 4/nc*ones(nc,1);
sigma = @(xy) 1; 
[pd2,h] = discrete_optimal_transport(disk,F,uv,sigma,area);

uv_new = compute_centroid(disk,pd2);
 