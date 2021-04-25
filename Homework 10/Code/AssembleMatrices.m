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

% Luca Wester <luca.wester@fau.de>, July 2019  

function [A,M] = AssembleMatrices(coord,elemNodeTable,coeff_a,coeff_r)

% Start of matrix assembly
n_vertices = size(coord, 1); 
n_elem = size(elemNodeTable, 1);

A  = sparse(n_vertices, n_vertices); %diffusion matrix
M  = sparse(n_vertices, n_vertices); %weighted mass matrix

% gradients of the basis functions in the reference element (constant!)
grd_bas_fcts = [ -1 -1 ; 1 0 ; 0 1 ]' ;
    
% We loop over the elements of the mesh,
% and add the contributions of each element to each matrix

% At each element we use the quadrature formula which uses 
% the function values at the midpoint of each side:
% \int_T  f  \approx  |T| ( f(a12) + f(a23) + f(a31) ) / 3.
% This formula is exact for quadratic polynomials

for el = 1 : n_elem

    %node indices of current element
    v_elem = elemNodeTable( el, : );
    
    %coordinates of the 3 corner nodes
    v1 = coord( v_elem(1), :)' ; 
    v2 = coord( v_elem(2), :)' ;
    v3 = coord( v_elem(3), :)' ; 
    
    %midpoints of the three sides
    a12 = (v1 + v2) / 2;
    a23 = (v2 + v3) / 2;
    a31 = (v3 + v1) / 2;
  
    % derivative of the affine transformation from the reference
    % element onto the current element
    D = [ v2-v1  v3-v1 ];
    
    % element area
    el_area = abs(det(D)) * 0.5;
    
    Dinv = inv(D);

    % quadrature rule applied to the diffusion coefficient
    coeff_a_quad = 1/3*(coeff_a(a12(1),a12(2)) + coeff_a(a23(1),a23(2)) + coeff_a(a31(1),a31(2)));
    
    % local contributions by the reactive part
    r12 = coeff_r(a12(1),a12(2));
    r23 = coeff_r(a23(1),a23(2));
    r31 = coeff_r(a31(1),a31(2));
    
    local_r = 1/3*[r12/4+r31/4, r12/4, r31/4; ...
                   r12/4,r12/4+r23/4,r23/4; ...
                   r31/4,r23/4,r31/4+r23/4];
    
    % computation of the element matrices
    el_mat_A = coeff_a_quad * grd_bas_fcts' * (Dinv*Dinv') * grd_bas_fcts * el_area;
    el_mat_M = el_area * local_r;
     
    % contributions added to the global matrices
    A( v_elem, v_elem ) = A( v_elem, v_elem ) + el_mat_A;
    M( v_elem, v_elem ) = M( v_elem, v_elem ) + el_mat_M;
    
end


end
