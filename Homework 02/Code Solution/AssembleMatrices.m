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

    for element_index = 1:n_elem
        vertex_indices = elemNodeTable(element_index,:);
        v = coord(vertex_indices,:); % 3x2. Vertices of the triangle
        B2C = [v'; 1,1,1]; % barycentric to cartesian coord. matrix
        C2B = inv(B2C);    % cartesian to barycentric coord. matrix
        area = 0.5*abs(det(B2C));
        phi = cell(3);
        for i = 1:3
            % In barycentric coordinates, phi_i is just lambda_i.
            % In cartesian coordinates, phi_i is:
            phi{i} = @(x1,x2) C2B(i,1)*x1 + C2B(i,2)*x2 + C2B(i,3);
        end
        grad_phi = C2B(:,1:2); % 3x2. Row i is the gradient of phi_i.
        for i = 1:3
            for j = 1:3
                % This code is optimized for readability, not speed.
                % Inlining the anonymous functions and the quadrature
                % routine would make everything at least 10x faster.
                integrandA = @(x1,x2) coeff_a(x1,x2)*dot(grad_phi(i,:),grad_phi(j,:));
                integrandM = @(x1,x2) coeff_r(x1,x2)*phi{i}(x1,x2)*phi{j}(x1,x2);
                integralA = trapezoidal_quadrature(v, area, integrandA);
                integralM = trapezoidal_quadrature(v, area, integrandM);
                vi_idx = vertex_indices(i);
                vj_idx = vertex_indices(j);
                A(vi_idx, vj_idx) = A(vi_idx, vj_idx) + integralA;
                M(vi_idx, vj_idx) = M(vi_idx, vj_idx) + integralM;
            end
        end
    end
end













