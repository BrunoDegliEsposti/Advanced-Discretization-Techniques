function [uh] = TimeStepping(coeff_a, coeff_r, u_ex, f, coord, elemNodeTable, boundary, deltat, T, theta)
    n_vertices = size(coord,1);
    n_boundary_edges = size(boundary,1);
    
    [A,wM] = AssembleMatrices(coord,elemNodeTable,coeff_a,coeff_r);
    M = unweightedMassMatrix(coord,elemNodeTable);
    A = A + wM;
    S = M + deltat*theta*A;
    W = M - deltat*(1-theta)*A;
    % With this notation, theta-method iterations will be:
    % S u_{k+1} = W u_k + deltat (theta rhs_{k+1} + (1-theta) rhs_k)
    
    for i = 1:n_boundary_edges
        if boundary(i,1) == 1
            vertex_id = boundary(i,2);
            S(vertex_id,:) = 0;
            S(vertex_id,vertex_id) = 1;
            vertex_id = boundary(i,3);
            S(vertex_id,:) = 0;
            S(vertex_id,vertex_id) = 1;
        end
    end
    
    t = 0;
    nsteps = round(T/deltat);
    uh = zeros(n_vertices, nsteps+1);
    rhs = zeros(n_vertices, nsteps+1);
    
    uh(:,1) = u_ex(coord(:,1),coord(:,2),0);
    for k = 1:nsteps+1
        rhs(:,k) = AssembleRHS(coord,elemNodeTable,f,deltat*(k-1));
    end
    
    for k = 1:nsteps
        r = W*uh(:,k) + deltat * (theta*rhs(:,k+1) + (1-theta)*rhs(:,k));
        for i = 1:n_boundary_edges
            if boundary(i,1) == 1
                vertex_id = boundary(i,2);
                vertex = coord(vertex_id,:);
                r(vertex_id) = u_ex(vertex(1),vertex(2),t+deltat);
                vertex_id = boundary(i,3);
                vertex = coord(vertex_id,:);
                r(vertex_id) = u_ex(vertex(1),vertex(2),t+deltat);
            end
        end
        uh(:,k+1) = S\r;
        t = t + deltat;
    end
end

