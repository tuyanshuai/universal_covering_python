 
addpath('mfile')

%% Example 1: Face

[F, V]= read_off('../data/faceg0.off');

map = disk_conformal_map(V,F);

plot_mesh(map,F); view([-90 90]);

% evaluate the angle distortion
angle_distortion(V,F,map);
 
% 
% %%
% 
% load('human_brain.mat')
% 
% % plot_mesh(v,f); view([90 0]);
% % can also include the third input if an additional quantity is defined on vertices
% plot_mesh(v,f,mean_curv); view([90 0]); 
% 
% map = disk_conformal_map(v,f);
% 
% % plot_mesh(map,f); 
% % can also include the third input if an additional quantity is defined on vertices
% plot_mesh(map,f,mean_curv); 
% 
% % evaluate the angle distortion
% angle_distortion(v,f,map);