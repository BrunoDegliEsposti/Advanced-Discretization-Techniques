%Assembles the RHS of a Lagrange-P1 finite element system given a
%forcing f

%Expected input:
%coord - nc x 2 matrix, storing the coordinates of each node
%elemNodeTable - nt x 3 matrix, storing the element connectivity information
%f = @(x,y) ... - function handle, forcing term 
%t - point in time
 
% Luca Wester <luca.wester@fau.de>, November 2019 

function rhs = AssembleRHS(coord,elemNodeTable,f)

% Start of the right-hand side assembly
n_vertices = size(coord, 1); 
n_elem = size(elemNodeTable, 1);

rhs = zeros(n_vertices, 1);

%Linear or quadratic elements?
if(size(elemNodeTable,2) == 3)
  el_type = 1;
elseif(size(elemNodeTable,2) == 6)
  el_type = 2;
end
    
% We loop over the elements of the mesh,
% and add the contributions of each element to the right-hand side rhs

% At each element we use the quadrature formula which uses 
% the function values at the midpoint of each side:
% \int_T  f  \approx  |T| ( f(a12) + f(a23) + f(a31) ) / 3.
% This formula is exact for quadratic polynomials

%Linear Lagrange (no different from Programming exercise 1)
if(el_type == 1)

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
      f12 = f(a12(1),a12(2));  f23 = f(a23(1),a23(2));  f31 = f(a31(1),a31(2));
    
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

%Quadratic Lagrange
elseif(el_type == 2)

  for el = 1 : n_elem

      %node indices of current element
      v_elem = elemNodeTable( el, : );
      
      %coordinates of the 3 corner nodes
      v1 = coord( v_elem(1), :)' ; 
      v2 = coord( v_elem(2), :)' ;
      v3 = coord( v_elem(3), :)' ; 
      
      %coordinates of the 3 edge midpoints
      a23 = coord( v_elem(4), :)' ; 
      a31 = coord( v_elem(5), :)' ;
      a12 = coord( v_elem(6), :)' ; 
      
      % evaluation of f at the quadrature points
      f12 = f(a12(1),a12(2));  f23 = f(a23(1),a23(2));  f31 = f(a31(1),a31(2));
    
      % derivative of the affine transformation from the reference
      % element onto the current element
      D = [ v2-v1  v3-v1 ];
      
      % element area
      el_area = abs(det(D)) * 0.5;

      % computation of the element load vector
      % note: Like with the mass matrix, phi_i(x_j) != 0 if and only if i = j
      % edge midpoint quadrature rule thus leads to an element load vector like:
      f_el = [0;0;0;f23;f31;f12] * (el_area/3);
      
      % contributions added to the global load vector
      rhs( v_elem ) = rhs( v_elem ) + f_el;
      
  end

end

end
