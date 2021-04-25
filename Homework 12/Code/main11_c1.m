% PROGRAMMING EXERCISE 11: CONVECTION-DOMINATED PROBLEMS WITH SUPG
% Bruno Degli Esposti, Xingyu Xu
% 21/01/20 - 04/02/20
% Code tested in MATLAB only

% Comment on the results:
% As we decrease the value of delta_T, the oscillation artifacts
% in the solution become larger and larger, and in the limit
% as delta_T -> 0 we get the same result as the implementation
% without SUPG. On the other hand, if we increase the value
% of delta_T too much, we get significant artificial diffusion
% along the direction of the advection field c(x,y)
% (but not the other, as expected), and so the peak value of 1
% quickly spreads around and decays.
% A good choice for delta_T balances these two opposite
% and undesirable effects.
% Here are some values that gave us satisfactory results:
% dx = 0.05     dt = 0.1     delta_T = 0.3 * diameter(T)
% dx = 0.05     dt = 0.03    delta_T = 0.2 * diameter(T)
% dx = 0.025    dt = 0.1     delta_T = 0.5 * diameter(T)
% dx = 0.025    dt = 0.03    delta_T = 0.1 * diameter(T)
% dx = 0.01     dt = 0.1     delta_T = 1.0 * diameter(T)
% dx = 0.01     dt = 0.03    delta_T = 0.2 * diameter(T)

%% Parameters

% of the PDE
coeff_a = @(x,y) 0.0001;
coeff_c = @(x,y) [0.5; 0.5];
coeff_r = @(x,y) 0;
uD = @(x,y,t) 0;
u0 = @(x,y) sin(5*pi*x) .* sin(5*pi*y) .* (x>=0.2 & x<=0.4 & y>=0.2 & y<= 0.4);
f = @(x,y,t) 0;

% of the space discretization
dx = 0.025;
N = round(1/dx);
BC = [1;1;1;1];
[coord,elemNodeTable,boundary] = gen_mesh_rectangle(N,N,0,1,0,1,BC);
n_vertices = size(coord,1);
n_elem = size(elemNodeTable, 1);

% of the time discretization
dt = 0.03;
T = 1;
nsteps = round(T/dt);

% of the SUPG method
deltaT_SUPG = zeros(n_elem,1);
for i = 1:n_elem
    v_elem = elemNodeTable(i,:);
    v1 = coord(v_elem(1),:)';
    v2 = coord(v_elem(2),:)';
    v3 = coord(v_elem(3),:)';
    l1 = norm(v3-v2);
    l2 = norm(v1-v3);
    l3 = norm(v2-v1);
    deltaT_SUPG(i) = max([l1,l2,l3]);
end
deltaT_SUPG = 0.1*deltaT_SUPG;

%% Solution

uh = CrankNicolson(coeff_a, coeff_c, coeff_r, f, uD, u0, ...
    coord, elemNodeTable, boundary, dt, T);
uh_SUPG = CrankNicolson_SUPG(coeff_a, coeff_c, coeff_r, f, uD, u0, ...
    coord, elemNodeTable, boundary, deltaT_SUPG, dt, T);

%% Plots

figure(1);
for i = 1:3
    for k = 0:nsteps
        subplot(1,2,1);
        trisurf(elemNodeTable,coord(:,1),coord(:,2),uh(:,k+1));
        zlim([-0.1,1]); shading interp;
        
        subplot(1,2,2);
        trisurf(elemNodeTable,coord(:,1),coord(:,2),uh_SUPG(:,k+1));
        zlim([-0.1,1]); shading interp;
        
        drawnow; pause(0.01);
    end
end


