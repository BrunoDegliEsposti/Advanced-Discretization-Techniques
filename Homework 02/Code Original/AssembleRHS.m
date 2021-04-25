%Assembles the RHS of a Lagrange-P1 finite element system given a
%forcing f

%Expected input:
%coord - nc x 2 matrix, storing the coordinates of each node
%elemNodeTable - nt x 3 matrix, storing the element connectivity information
%f = @(x,y) ... - function handle, forcing term 

function rhs = AssembleRHS(coord,elemNodeTable,f)

% Start of the right-hand side assembly
n_vertices = size(coord, 1); 
n_elem = size(elemNodeTable, 1);

rhs = zeros(n_vertices, 1);

for el = 1 : n_elem

  %TODO: Fill with life
    
end


end
