%Tests AssembleMatrices.m in the case of the stationary elliptic problem

%  - nabla \cdot (a nabla u) + ru  = f      in Omega
%                                u = u_D    on Gamma

% ############################## Test case ################################
%u
u_ex = @(x,y) sin(2*pi*x) .* exp(y);

%\nabla u
grad_u_ex_1 = @(x,y) 2*pi*cos(2*pi*x).*exp(y);
grad_u_ex_2 = @(x,y) u_ex(x,y);

%\Delta u
lap_u = @(x,y) (-4*pi*pi + 1)*u_ex(x,y);

%a, r
coeff_a = @(x,y) x+sin(y).^2;
coeff_r = @(x,y) 2.*x+2.*y;

%\nabla a
d1a = @(x,y) 1;
d2a = @(x,y) 2*sin(y).*cos(y);

%f
f = @(x,y) -(d1a(x,y).*grad_u_ex_1(x,y) + d2a(x,y).*grad_u_ex_2(x,y) +  ...
           coeff_a(x,y).*lap_u(x,y)) + coeff_r(x,y).*u_ex(x,y);

% #########################################################################

for N = [4,8,16,32]
    [coord,elemNodeTable,boundary] = gen_mesh_rectangle(N, N, ...
        0, 1, 0, 2, [1;1;1;1], 2);
    [A,M] = AssembleMatrices(coord,elemNodeTable,coeff_a,coeff_r);
    B = A+M;
    rhs = AssembleRHS(coord,elemNodeTable,f);
    n_boundary_edges = size(boundary,1);
    % The essential boundary conditions are introduced by the
    % diagonalization technique (ref. section A.2.3 of 'Solving
    % Numerical PDEs' - Formaggia, Saleri, Veneziani).
    for i = 1:n_boundary_edges
        if boundary(i,1) == 1 % Dirichlet BC
            for j = 2:size(boundary,2)
                vertex_id = boundary(i,j);
                vertex = coord(vertex_id,:);
                B(vertex_id,:) = 0;
                B(vertex_id,vertex_id) = 1;
                rhs(vertex_id) = u_ex(vertex(1),vertex(2));
            end
        end
    end
    uh = B\rhs;
    u = u_ex(coord(:,1),coord(:,2));
    e = uh-u;
    % Check if the error goes to zero as N increases
    fprintf('N = %d, err max = %f\n',N,norm(e,Inf));
end

% We end with two plots
figure(1);
subplot(1,2,1), trisurf(elemNodeTable(:,1:3), coord(:,1), coord(:,2), uh), ...
    shading interp;
subplot(1,2,2), trisurf(elemNodeTable(:,1:3), coord(:,1), coord(:,2), u), ...
    shading interp;
figure(2);
trisurf(elemNodeTable(:,1:3),coord(:,1),coord(:,2),u-uh), shading interp;








