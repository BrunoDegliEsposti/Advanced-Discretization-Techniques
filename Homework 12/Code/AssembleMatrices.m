% PROGRAMMING EXERCISE 11: CONVECTION-DOMINATED PROBLEMS WITH SUPG
% Bruno Degli Esposti, Xingyu Xu
% 21/01/20 - 04/02/20
% Code tested in MATLAB only

function [A,B,R] = AssembleMatrices(coord,elemNodeTable,...
        coeff_a,coeff_c,coeff_r)
    %ASSEMBLEMATRICES Assemble diffusion, advection, reaction
    % matrices using Lagrange P1 elements
    n_vertices = size(coord, 1);
    n_elem = size(elemNodeTable, 1);
    
    A = sparse(n_vertices, n_vertices); % diffusion matrix
    B = sparse(n_vertices, n_vertices); % advection matrix
    R = sparse(n_vertices, n_vertices); % reaction/w.mass matrix
    
    % basis functions N_1 = 1-x-y, N_2 = x, N_3 = y evaluated
    % in the quadrature nodes a4 = (0.5,0.5), a5 = (0,0.5), a6 = (0.5,0)
    N = [0.0, 0.5, 0.5;
         0.5, 0.0, 0.5;
         0.5, 0.5, 0.0];
    % gradients of the basis functions N_i in the reference element
    DN = [-1, 1, 0;
          -1, 0, 1];
    
    for i = 1:n_elem
        v_elem = elemNodeTable(i,:);
        % coordinates of the 3 corner nodes
        v1 = coord(v_elem(1),:)';
        v2 = coord(v_elem(2),:)';
        v3 = coord(v_elem(3),:)';
        % midpoints of the three sides
        a4 = 0.5 * (v2+v3);
        a5 = 0.5 * (v3+v1);
        a6 = 0.5 * (v1+v2);
        % derivative of the affine transformation from the reference
        % element onto the current element
        D = [v2-v1, v3-v1];
        Dinv = inv(D);
        el_area = 0.5 * abs(det(D));
        
        % evaluate the coefficients a(x,y), c(x,y), r(x,y)
        % in the quadrature nodes a4, a5, a6:
        % a is a scalar
        a = coeff_a(a6(1),a6(2)) + ...
            coeff_a(a4(1),a4(2)) + ...
            coeff_a(a5(1),a5(2));
        % c is a 2x3 matrix
        c = [coeff_c(a4(1),a4(2)), ...
             coeff_c(a5(1),a5(2)), ...
             coeff_c(a6(1),a6(2))];
        % r is a 3x3 matrix
        r = diag([coeff_r(a4(1),a4(2)),...
                  coeff_r(a5(1),a5(2)),...
                  coeff_r(a6(1),a6(2))]);
        
        % vectorized computation of the 3x3 local matrix
        A_local = (1/3) * el_area * a * DN' * (Dinv*Dinv') * DN;
        B_local = (1/3) * el_area * N * c' * Dinv' * DN;
        R_local = (1/3) * el_area * N' * r * N;
        
        % contributions added to the global matrices
        A(v_elem, v_elem) = A(v_elem, v_elem) + A_local;
        B(v_elem, v_elem) = B(v_elem, v_elem) + B_local;
        R(v_elem, v_elem) = R(v_elem, v_elem) + R_local;
    end
end













