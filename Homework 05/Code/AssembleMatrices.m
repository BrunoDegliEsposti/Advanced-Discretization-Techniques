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
    n_vertices = size(coord, 1);
    n_elem = size(elemNodeTable, 1);
    A  = sparse(n_vertices, n_vertices);
    M  = sparse(n_vertices, n_vertices);
    if size(elemNodeTable,2) == 3 % Linear Lagrange elements
        % Start of matrix assembly
        grd_bas_fcts = [ -1 -1 ; 1 0 ; 0 1 ]' ;        
        for el = 1:n_elem
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
    elseif size(elemNodeTable,2) == 6 % Quadratic Lagrange elements
        for el = 1:n_elem
            v_elem = elemNodeTable(el,:); % 6 vertices a1-a6
            a1 = coord(v_elem(1),:)';
            a2 = coord(v_elem(2),:)';
            a3 = coord(v_elem(3),:)';
            a4 = coord(v_elem(4),:)'; % midpoint of a2 and a3
            a5 = coord(v_elem(5),:)'; % midpoint of a3 and a1
            a6 = coord(v_elem(6),:)'; % midpoint of a1 and a2
            a_a4 = coeff_a(a4(1),a4(2));
            a_a5 = coeff_a(a5(1),a5(2));
            a_a6 = coeff_a(a6(1),a6(2));
            r_a4 = coeff_r(a4(1),a4(2));
            r_a5 = coeff_r(a5(1),a5(2));
            r_a6 = coeff_r(a6(1),a6(2));
            D = [a2-a1, a3-a1];
            el_area = 0.5*abs(det(D));
            Dinv = inv(D);
            DinvDTinv = Dinv*Dinv';
            N = cell(6);
            N{1} = @(x,y) (1-x-y)*(1-2*x-2*y);
            N{2} = @(x,y) x*(2*x-1);
            N{3} = @(x,y) y*(2*y-1);
            N{4} = @(x,y) 4*x*y;
            N{5} = @(x,y) 4*y*(1-x-y);
            N{6} = @(x,y) 4*x*(1-x-y);
            gradN = cell(6);
            gradN{1} = @(x,y) [4*x+4*y-3; 4*x+4*y-3];
            gradN{2} = @(x,y) [4*x-1; 0];
            gradN{3} = @(x,y) [0; 4*y-1];
            gradN{4} = @(x,y) [4*y; 4*x];
            gradN{5} = @(x,y) [-4*y; -4*x-8*y+4];
            gradN{6} = @(x,y) [-8*x-4*y+4; -4*x];
            el_mat_A = zeros(6,6);
            el_mat_M = zeros(6,6);
            for i = 1:6
                for j = 1:i
                    el_mat_A(i,j) = (1/3)*el_area* ...
                        (a_a4*dot(gradN{i}(0.5,0.5),DinvDTinv*gradN{j}(0.5,0.5)) ...
                        +a_a5*dot(gradN{i}(0.0,0.5),DinvDTinv*gradN{j}(0.0,0.5)) ...
                        +a_a6*dot(gradN{i}(0.5,0.0),DinvDTinv*gradN{j}(0.5,0.0)));
                    el_mat_M(i,j) = (1/3)*el_area* ...
                        (r_a4*N{i}(0.5,0.5)*N{j}(0.5,0.5) ...
                        +r_a5*N{i}(0.0,0.5)*N{j}(0.0,0.5) ...
                        +r_a6*N{i}(0.5,0.0)*N{j}(0.5,0.0));
                end
            end
            % exploit symmetry
            el_mat_A = el_mat_A + tril(el_mat_A,-1)';
            el_mat_M = el_mat_M + tril(el_mat_M,-1)';
            A(v_elem, v_elem) = A(v_elem, v_elem) + el_mat_A;
            M(v_elem, v_elem) = M(v_elem, v_elem) + el_mat_M;
        end
    else
        error('Cannot assemble matrices: mesh not supported')
    end
end












