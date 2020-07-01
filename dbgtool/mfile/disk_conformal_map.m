function map = disk_conformal_map(vertex,face)

% Compute the disk conformal mapping using the method in [1].
%
% Input:
% vertex: nv x 3 vertex coordinates of a simply-connected open triangle mesh
% face: nf x 3 triangulations of a simply-connected open triangle mesh
% 
% Output:
% map: nv x 2 vertex coordinates of the disk conformal parameterization
% 
% Remark:
% 1. Please make sure that the input mesh does not contain any 
%    unreferenced vertices/non-manifold vertices/non-manifold edges.
% 2. Please remove all valence 1 boundary vertices (i.e. vertices with 
%    only 1 face attached to them) before running the program.
% 
% If you use this code in your own work, please cite the following paper:
% [1] P. T. Choi and L. M. Lui, 
%     "Fast Disk Conformal Parameterization of Simply-Connected Open Surfaces."
%     Journal of Scientific Computing, 65(3), pp. 1065-1090, 2015.
 
parameter.north = 5;
parameter.south = 100;
parameter.threshold = 0.00001;

bdy_index = compute_bd(face);
map = disk_harmonic_map(face,vertex);
z =  map(:,1) + 1i*map(:,2);
fprintf('Initialization completed.\n')
%% North Pole iteration
% Use the Cayley transform to map the disk to the upper half plane
% All boundary points will be mapped to the real line

mu = beltrami_coefficient(map, face, vertex); 
mu_v = f2v(vertex,face)*mu;
bdy_index_temp = [bdy_index(2:end);bdy_index(1)];
[~, least] = min(abs(mu_v(bdy_index))+abs(mu_v(bdy_index_temp)));
z = z* exp(-1i*(angle(z(bdy_index(least)))+angle(z(bdy_index(mod(least,length(bdy_index))+1)))/2));
g = 1i*(1 + z)./(1 - z);

% fix the points near the puncture, i.e. near z = 1
[~, ind] = sort(-real(z));
fixed = setdiff(ind(1:max(round(length(vertex)/parameter.north),min(100,length(z)))), bdy_index);
fixed = [fixed; find(real(g) == max(real(g))); find(real(g) == min(real(g)))];

P = [real(g),imag(g),ones(length(g),1)];
mu = beltrami_coefficient(P, face, vertex); 

% compute the updated x coordinates
target = P(fixed,1);
A = generalized_laplacian(P,face,mu); Ax = A; Ay = A;
b = -Ax(:,fixed)*target;
b(fixed) = target;
Ax(fixed,:) = 0; Ax(:,fixed) = 0;
Ax = Ax + sparse(fixed,fixed,ones(length(fixed),1), size(A,1), size(A,2));
x = Ax\b;

% compute the updated y coordinates
target = P(fixed,2); 
fixed = [fixed;bdy_index];
target = [target;zeros(length(bdy_index),1)];
b = -Ay(:,fixed)*target;
b(fixed) = target;
Ay(fixed,:) = 0; Ay(:,fixed) = 0;
Ay = Ay + sparse(fixed,fixed,ones(length(fixed),1), size(Ay,1), size(A,2));
y = Ay\b;

g_new = complex(x,y);
z_new = (g_new - 1i)./(g_new + 1i);
map = [real(z_new), imag(z_new)];

if sum(sum(isnan(map))) ~= 0 
    % use the old result in case of getting NaN entries
    map = map_prev;
    z_new = map(:,1) + 1i*map(:,2);
end

fprintf('North pole step completed.\n')

%% Reflection along the unit circle

f_temp = face + length(vertex);
a = sort(bdy_index + length(vertex));
for i = length(a):-1:1
    f_temp(f_temp == a(i)) = a(i) - length(vertex);
    f_temp = f_temp - (f_temp >a(i));
end
f_filled = [face;fliplr(f_temp)];

z_filled = [z_new;1./conj(z_new)];
z_filled(bdy_index + length(vertex)) = [];

energy_old = 0;
energy = mean(abs(beltrami_coefficient([real(z_new),imag(z_new),0.*z_new], face, vertex)));

iteration_count = 1; 

map_opt = map;

fprintf('Reflection completed.\n')

%% South pole iteration
% Iteratively compose the map with a quasi-conformal map,
% then normalize the boundary
while abs(energy_old-energy) > parameter.threshold

    energy_old = energy;
    
    mu = beltrami_coefficient([real(z_new),imag(z_new),0.*z_new], face, vertex);
    mu_filled = [mu;1/3*((z_new(face(:,1))./(conj(z_new(face(:,1))))).^2 + ...
        (z_new(face(:,2))./(conj(z_new(face(:,2))))).^2 + ...
        (z_new(face(:,3))./(conj(z_new(face(:,3))))).^2).*conj(mu)./...
        abs(((z_new(face(:,1))./(conj(z_new(face(:,1))))).^2 + ...
        (z_new(face(:,2))./(conj(z_new(face(:,2))))).^2 + ...
        (z_new(face(:,3))./(conj(z_new(face(:,3))))).^2))];

    % fix the points near infinity
    [~, ind] = sort(-abs(z_filled));
    fixed2 = ind(1:max(round(length(z_filled)/parameter.south),...
        min(100,length(z))));
    map_filled = linear_beltrami_solver(...
        [real(z_filled),imag(z_filled),0.*z_filled],f_filled,mu_filled,...
        fixed2,[real(z_filled(fixed2)), imag(z_filled(fixed2))]);

    z_big = complex(map_filled(:,1),map_filled(:,2));
    z_final = z_big(1:length(vertex));
    
    % normalization
    z_final = z_final - mean(z_final); % move centroid to zero
    if max(abs(z_final))>1
        z_final = z_final/(max(abs(z_final))); % map it into unit circle
    end
    mu_temp = beltrami_coefficient([real(z_final),imag(z_final),0.*z_final],face,vertex);
    map_temp = linear_beltrami_solver(...
        [real(z_final),imag(z_final),0.*z_final],face,mu_temp,...
        bdy_index,[real(z_final(bdy_index)./abs(z_final(bdy_index))), ...
        imag(z_final(bdy_index)./abs(z_final(bdy_index)))]);

    z_new = map_temp(:,1) + 1i*map_temp(:,2);

    z_filled = [z_new;1./conj(z_new)];
    z_filled(bdy_index + length(vertex)) = [];

    map = [real(z_new), imag(z_new)];
    
    if sum(sum(isnan(map))) ~= 0 
        % use the previous result in case of getting NaN entries
        map = map_opt;
        map(:,1) = -map(:,1);
        return;
    end

    energy = mean(abs(beltrami_coefficient(map, face, vertex)));
    map_opt = map;
  
    fprintf('Iteration %d: mean(|mu|) = %.4f\n',[iteration_count,energy]);
    iteration_count = iteration_count+1;
    
    if iteration_count > 5
        % it usually converges within 5 iterations so we set 5 here
        break;
    end
end
map = map_opt;
map(:,1) = -map(:,1);

fprintf('South pole step completed.\n')