%% PROGRAMMING EXERCISE 5: STOKES WITH TAYLOR-HOOD ELEMENTS
% Bruno Degli Esposti, Xingyu Xu
% 19/11/19 - 3/12/19
% Code tested in MATLAB only

%% Parameters

n = 30;         % horizontal mesh subdivisions
m = 30;         % vertical mesh subdivisions
g = 9.81;       % gravity constant
umax = 0.5;     % horizontal vel. of the fluid at the top of the cavity

%% Assembly

[coord,elemNodeTable,boundary] = gen_mesh_rectangle(n, m, ...
    0, 1, 0, 1, [1;1;1;1], 2);
[A,~] = AssembleMatrices(coord, elemNodeTable, @(x,y) 1, @(x,y) 1);
[B1,B2] = AssembleMixedMatrix(coord, elemNodeTable);
[N1,N2] = size(B1);
S = [A, zeros(N2), B1'; zeros(N2), A, B2'; B1, B2, zeros(N1)];
g_rhs = AssembleRHS(coord, elemNodeTable, @(x,y) g);
rhs = [zeros(N2,1); -g_rhs; zeros(N1,1)];

%% Boundary conditions

n_boundary_edges = size(boundary,1);
% The essential boundary conditions are introduced by the
% diagonalization technique (ref. section A.2.3 of 'Solving
% Numerical PDEs' - Formaggia, Saleri, Veneziani).
for i = 1:n_boundary_edges
    if boundary(i,1) == 1 % Dirichlet BC
        for j = 2:4
            vertex_id = boundary(i,j);
            vertex = coord(vertex_id,:);
            S(vertex_id,:) = 0;
            S(vertex_id+N2,:) = 0;
            S(vertex_id,vertex_id) = 1;
            S(vertex_id+N2,vertex_id+N2) = 1;
            if abs(vertex(2)-1) < 1e-14
                rhs(vertex_id) = umax;
                rhs(vertex_id+N2) = 0;
            else
                rhs(vertex_id) = 0;
                rhs(vertex_id+N2) = 0;
            end
        end
    end
end
% set p to 0 on the bottom left corner of the cavity
S(2*N2+1,:) = 0;
S(2*N2+1,2*N2+1) = 1;
rhs(2*N2+1) = 0;

%% Solution of the linear system and plots

sol = S\rhs;
u1 = sol(1:N2);
u2 = sol(N2+1:2*N2);
p = sol(2*N2+1:end);

figure(1);
hold on;
plot([0,0,1,1],[1,0,0,1],'k','LineWidth',4);
contourf(linspace(0,1,n+1), linspace(0,1,m+1), reshape(p,n+1,m+1)', 20);
hold off;

figure(2);
hold on;
plot([0,0,1,1],[1,0,0,1],'k','LineWidth',2);
quiver(coord(:,1), coord(:,2), u1, u2, 3);
hold off;

%% Comments on the results

% if umax = 0, the fluid is at rest and p is determined by the law
% of hydrostatic pressure p = p0 + rho*g*h, so the isolines are horizontal.

% As umax increses, the pressure start to vary in the horizontal
% direction as well, and this variation is greater at the top of the cavity.
% The pressure in the top corners (0,1) and (1,1) is singular:
% as n,m -> inf, the pressure goes to -inf and +inf (respectively).

% The velocity field has nonzero curl, and rotates around a point
% near (0.5,0.75). The center of rotation doesn't seem to be affected
% by the choice of parameters g and umax.





