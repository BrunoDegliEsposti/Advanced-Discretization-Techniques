%Assembles the RHS of a Lagrange-P1 finite element system given a
%forcing f at time t

%Expected input:
%coord - nc x 2 matrix, storing the coordinates of each node
%elemNodeTable - nt x 3 matrix, storing the element connectivity information
%f = @(x,y,t) ... - function handle, forcing term
%t - time at which f is evaluated
 
% Luca Wester <luca.wester@fau.de>, July 2019 

function rhs = AssembleRHSt(coord,elemNodeTable,f,t)

% Start of the right-hand side assembly
n_vertices = size(coord, 1); 
n_elem = size(elemNodeTable, 1);

rhs = zeros(n_vertices, 1);
    
% We loop over the elements of the mesh,
% and add the contributions of each element to the right-hand side rhs

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

    % evaluation of f at the quadrature points
    f12 = f(a12(1),a12(2),t);  f23 = f(a23(1),a23(2),t);  f31 = f(a31(1),a31(2),t);
  
    % derivative of the affine transformation from the reference
    % element onto the current element
    D = [ v2-v1  v3-v1 ];
    
    % element area
    el_area = abs(det(D)) * 0.5;

    % computation of the element load vector
    % note: basis function equal to 0.5 on two adjacent edges, 0 on the opposite edge
    f_el = [ (f12+f31)*0.5 ; (f12+f23)*0.5 ; (f23+f31)*0.5 ] * (el_area/3);
    
    % contributions added to the global load vector
    rhs( v_elem ) = rhs( v_elem ) + f_el;
    
end


end
