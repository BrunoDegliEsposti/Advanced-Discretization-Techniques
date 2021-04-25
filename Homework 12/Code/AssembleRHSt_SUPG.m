% PROGRAMMING EXERCISE 11: CONVECTION-DOMINATED PROBLEMS WITH SUPG
% Bruno Degli Esposti, Xingyu Xu
% 21/01/20 - 04/02/20
% Code tested in MATLAB only

function F_SUPG = AssembleRHSt_SUPG(coord, elemNodeTable, ...
        deltaT_SUPG, coeff_c, f, t)
    %ASSEMBLERHSTSUPG Assembles the SUPG term for the RHS
    % of a Lagrange-P1 finite element system given a forcing f at time t
    n_vertices = size(coord, 1);
    n_elem = size(elemNodeTable, 1);
    
    F_SUPG = zeros(n_vertices, 1);
    
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
        
        % evaluate the coefficient c(x,y) and the forcing f(x,y,t)
        % in the quadrature nodes a4, a5, a6 at time t:
        % c is a 2x3 matrix
        c = [coeff_c(a4(1),a4(2)), ...
             coeff_c(a5(1),a5(2)), ...
             coeff_c(a6(1),a6(2))];
        % forcing is a 3x1 vector
        forcing = [f(a4(1),a4(2),t);
                   f(a5(1),a5(2),t);
                   f(a6(1),a6(2),t)];
        
        % vectorized computation of the local contribution to F_SUPG
        F_SUPG_local = (1/3) * el_area * DN' * Dinv * c * forcing;
        
        % contribution added to the global matrix
        F_SUPG(v_elem) = F_SUPG(v_elem) + deltaT_SUPG(i)*F_SUPG_local;
    end
end
