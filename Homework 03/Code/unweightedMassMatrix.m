function [M] = unweightedMassMatrix(coord,elemNodeTable)
    % Start of matrix assembly
    n_vertices = size(coord, 1); 
    n_elem = size(elemNodeTable, 1);

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
        for i = 1:3
            for j = 1:3
                integrandM = @(x1,x2) phi{i}(x1,x2)*phi{j}(x1,x2);
                integralM = trapezoidal_quadrature(v, area, integrandM);
                vi_idx = vertex_indices(i);
                vj_idx = vertex_indices(j);
                M(vi_idx, vj_idx) = M(vi_idx, vj_idx) + integralM;
            end
        end
    end
end
