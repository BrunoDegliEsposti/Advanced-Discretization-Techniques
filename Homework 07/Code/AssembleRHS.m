%Assembles the RHS of a Lagrange-P1 finite element system given a
%forcing f

%Expected input:
%coord - nc x 2 matrix, storing the coordinates of each node
%elemNodeTable - nt x 3 matrix, storing the element connectivity information
%f = @(x,y) ... - function handle, forcing term

function rhs = AssembleRHS(coord,elemNodeTable,f)
    n_vertices = size(coord, 1);
    n_elem = size(elemNodeTable, 1);
    rhs = zeros(n_vertices, 1);
    
    if size(elemNodeTable,2) == 3 % Linear Lagrange elements
        % Start of the right-hand side assembly
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
    elseif size(elemNodeTable,2) == 6 % Quadratic Lagrange elements
        for el = 1:n_elem
            v_elem = elemNodeTable(el,:); % 6 vertices a1-a6
            a1 = coord(v_elem(1),:)';
            a2 = coord(v_elem(2),:)';
            a3 = coord(v_elem(3),:)';
            a4 = coord(v_elem(4),:)'; % midpoint of a2 and a3
            a5 = coord(v_elem(5),:)'; % midpoint of a3 and a1
            a6 = coord(v_elem(6),:)'; % midpoint of a1 and a2
            f_a4 = f(a4(1),a4(2));
            f_a5 = f(a5(1),a5(2));
            f_a6 = f(a6(1),a6(2));
            D = [a2-a1, a3-a1];
            el_area = 0.5*abs(det(D));
            N = cell(6);
            N{1} = @(x,y) (1-x-y)*(1-2*x-2*y);
            N{2} = @(x,y) x*(2*x-1);
            N{3} = @(x,y) y*(2*y-1);
            N{4} = @(x,y) 4*x*y;
            N{5} = @(x,y) 4*y*(1-x-y);
            N{6} = @(x,y) 4*x*(1-x-y);
            f_el = zeros(6,1);
            for i = 1:6
                f_el(i) = (1/3)*el_area*(f_a4*N{i}(0.5,0.5) ...
                                        +f_a5*N{i}(0.0,0.5) ...
                                        +f_a6*N{i}(0.5,0.0));
            end
            rhs(v_elem) = rhs(v_elem) + f_el;
        end
    else
        error('Cannot assemble RHS: mesh not supported')
    end
end












