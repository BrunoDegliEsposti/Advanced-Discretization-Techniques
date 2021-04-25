% PROGRAMMING EXERCISE 5
% Bruno Degli Esposti, Xingyu Xu
% 19/11/19 - 3/12/19
% Code tested in MATLAB only

function [B1,B2] = AssembleMixedMatrix(coord, elemNodeTable)
    N1 = max(max(elemNodeTable(:,1:3)));
    N2 = size(coord, 1);
    B1 = sparse(N1,N2);
    B2 = sparse(N1,N2);
    n_elem = size(elemNodeTable, 1);
    assert(size(elemNodeTable,2)==6);
    for el = 1:n_elem
        v_elem = elemNodeTable(el,:); % indices of the 6 vertices a1-a6
        a1 = coord(v_elem(1),:)';
        a2 = coord(v_elem(2),:)';
        a3 = coord(v_elem(3),:)';
        D = [a2-a1, a3-a1];
        el_area = 0.5*abs(det(D));
        DinvT = inv(D)';
        psi = cell(3);
        psi{1} = @(x,y) 1-x-y;
        psi{2} = @(x,y) x;
        psi{3} = @(x,y) y;
        dxphi = cell(6);
        dxphi{1} = @(x,y) 4*x+4*y-3;
        dxphi{2} = @(x,y) 4*x-1;
        dxphi{3} = @(x,y) 0;
        dxphi{4} = @(x,y) 4*y;
        dxphi{5} = @(x,y) -4*y;
        dxphi{6} = @(x,y) -8*x-4*y+4;
        dyphi = cell(6);
        dyphi{1} = @(x,y) 4*x+4*y-3;
        dyphi{2} = @(x,y) 0;
        dyphi{3} = @(x,y) 4*y-1;
        dyphi{4} = @(x,y) 4*x;
        dyphi{5} = @(x,y) -4*x-8*y+4;
        dyphi{6} = @(x,y) -4*x;
        el_mat_B1 = zeros(3,6);
        el_mat_B2 = zeros(3,6);
        for i = 1:3
            for j = 1:6
                el_mat_B1(i,j) = -(1/3)*el_area*...
                    ((DinvT(1,1)*dxphi{j}(0.5,0.5)+DinvT(1,2)*dyphi{j}(0.5,0.5))*psi{i}(0.5,0.5) ...
                    +(DinvT(1,1)*dxphi{j}(0.0,0.5)+DinvT(1,2)*dyphi{j}(0.0,0.5))*psi{i}(0.0,0.5) ...
                    +(DinvT(1,1)*dxphi{j}(0.5,0.0)+DinvT(1,2)*dyphi{j}(0.5,0.0))*psi{i}(0.5,0.0));
                el_mat_B2(i,j) = -(1/3)*el_area*...
                    ((DinvT(2,1)*dxphi{j}(0.5,0.5)+DinvT(2,2)*dyphi{j}(0.5,0.5))*psi{i}(0.5,0.5) ...
                    +(DinvT(2,1)*dxphi{j}(0.0,0.5)+DinvT(2,2)*dyphi{j}(0.0,0.5))*psi{i}(0.0,0.5) ...
                    +(DinvT(2,1)*dxphi{j}(0.5,0.0)+DinvT(2,2)*dyphi{j}(0.5,0.0))*psi{i}(0.5,0.0));
            end
        end
        B1(v_elem(1:3), v_elem) = B1(v_elem(1:3), v_elem) + el_mat_B1;
        B2(v_elem(1:3), v_elem) = B2(v_elem(1:3), v_elem) + el_mat_B2;
    end
end




