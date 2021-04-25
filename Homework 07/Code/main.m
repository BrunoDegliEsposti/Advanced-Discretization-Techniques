% PROGRAMMING EXERCISE 6: UZAWA’S ALGORITHM
% Bruno Degli Esposti, Xingyu Xu
% 3/12/19 - 10/12/19
% Code tested in MATLAB only

%% Parameters

n = 20;         % horizontal mesh subdivisions
m = 20;         % vertical mesh subdivisions
g = 9.81;       % gravity constant
umax = 0.5;     % horizontal vel. of the fluid at the top of the cavity

omega = 1.8;    % convergence parameter in Uzawa's algorithm
tol = 1e-4;     % absolute tollerance in Uzawa's algorithm (stopping cr.)
kmax = 100;     % maximum number of iterations in Uzawa's algorithm

%% Assembly

[coord,elemNodeTable,boundary] = gen_mesh_rectangle(n, m, ...
    0, 1, 0, 1, [1;1;2;1], 2);
[A,~] = AssembleMatrices(coord, elemNodeTable, @(x,y) 1, @(x,y) 1);
[B1,B2] = AssembleMixedMatrix(coord, elemNodeTable);
[N1,N2] = size(B1);
g_rhs = AssembleRHS(coord, elemNodeTable, @(x,y) g);

A = [A, zeros(N2); zeros(N2), A];
B = [B1,B2];
[~,M] = AssembleMatrices(coord(1:(n+1)*(m+1),:), elemNodeTable(:,1:3), @(x,y) 1, @(x,y) 1);
rhs = [zeros(N2,1); -g_rhs; zeros(N1,1)];

%% Boundary conditions

no_slip_nodes = unique(boundary(boundary(:,1)==1,2:4));
lid_nodes = unique(boundary(boundary(:,1)==2,2:4));
no_slip_nodes = setdiff(no_slip_nodes, lid_nodes);

for i = no_slip_nodes'
    % set horizontal velocity to 0
    A(:,i) = 0;
    B(:,i) = 0;
    A(i,:) = 0;
    A(i,i) = 1;
    rhs(i) = 0;
    % set vertical velocity to 0
    A(:,N2+i) = 0;
    B(:,N2+i) = 0;
    A(N2+i,:) = 0;
    A(N2+i,N2+i) = 1;
    rhs(N2+i) = 0;
end

for i = lid_nodes'
    % set horizontal velocity to umax
    rhs = rhs - umax * [A(:,i); B(:,i)];
    A(:,i) = 0;
    B(:,i) = 0;
    A(i,:) = 0;
    A(i,i) = 1;
    rhs(i) = umax;
    % set vertical velocity to 0
    A(:,N2+i) = 0;
    B(:,N2+i) = 0;
    A(N2+i,:) = 0;
    A(N2+i,N2+i) = 1;
    rhs(N2+i) = 0;
end

%% Solution of the linear system and plots

p0 = zeros(N1,1);
[u,p,k] = uzawa(A,B,M,rhs(1:2*N2,1),rhs(2*N2+1:end,1),p0,omega,tol,kmax);
u1 = u(1:N2);
u2 = u(N2+1:2*N2);

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

% In line 6, p^{k+1}-p^k is implicitly seen as an element of the
% dual space of X_h. We need to include M because we have identified X_h
% with its dual X_h' through the L^2 scalar product,
% not through the Euclidean one that we get when we fix the psi_i basis.
% Then <p^{k+1}-p^k, \cdot>_{L^2} = <M(p^{k+1}-p^k), \cdot>_{R^N1},
% because M_ij = <psi_j,psi_i>_{L^2}.

% Testing with g = 9.81, umax = 0.5, tol = 1e-4, we see that
% the size of the mesh has no effect whatsoever on the number
% of iterations k required to meet the tolerance.
% On the other hand, a good choice for omega (1.8) can almost halve
% the number of iterations compared to omega=1:
% omega    k
%  1.00   54
%  1.25   44
%  1.50   37
%  1.75   32
%  1.80   31
%  1.85   32
%  1.90   46
%  1.95   85




