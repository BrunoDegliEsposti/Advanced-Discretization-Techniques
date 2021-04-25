%Assembles the RHS of a Lagrange-P1 finite element system given a
%forcing f

%Expected input:
%coord - nc x 2 matrix, storing the coordinates of each node
%elemNodeTable - nt x 3 matrix, storing the element connectivity information
%f = @(x,y) ... - function handle, forcing term 

function rhs = AssembleRHS(coord,elemNodeTable,f,t)
    % Start of the right-hand side assembly
    n_vertices = size(coord, 1); 
    n_elem = size(elemNodeTable, 1);

    rhs = zeros(n_vertices, 1);

    for element_index = 1:n_elem
        vertex_indices = elemNodeTable(element_index,:);
        v = coord(vertex_indices,:); % 3x2. Vertices of the triangle
        B2C = [v'; 1,1,1]; % barycentric to cartesian coord. matrix
        C2B = inv(B2C);    % cartesian to barycentric coord. matrix
        area = 0.5*abs(det(B2C));
        for i = 1:3
            phi_i = @(x1,x2) C2B(i,1)*x1 + C2B(i,2)*x2 + C2B(i,3);
            integrand = @(x1,x2) f(x1,x2,t)*phi_i(x1,x2);
            integral = trapezoidal_quadrature(v, area, integrand);
            vi_idx = vertex_indices(i);
            rhs(vi_idx) = rhs(vi_idx) + integral;
        end
    end
end













