% Assembles the diffusion system matrix A as well as the weighted mass 
% matrix M corresponding to the elliptic system 

%  - nabla \cdot (a(x) nabla u) + r(x) u  = f(x)      in Omega
%                                       u = g_D       on Gamma

% with scalar-valued coefficients a and r

%Expected input:
%coord - nc x 2 matrix, storing the coordinates of each node
%elemNodeTable - nt x 3 matrix, storing the element connectivity information
%coeff_a = @(x,y) ... , function handle, diffusion coefficient
%coeff_r = @(x,y) ... , function handle, reaction coefficient

function [A,M] = AssembleMatrices(coord,elemNodeTable,coeff_a,coeff_r)

% Start of matrix assembly
n_vertices = size(coord, 1); 
n_elem = size(elemNodeTable, 1);

A  = sparse(n_vertices, n_vertices); %diffusion matrix
M  = sparse(n_vertices, n_vertices); %weighted mass matrix


for el = 1 : n_elem

   %TODO: Fill with life
    
end


end
