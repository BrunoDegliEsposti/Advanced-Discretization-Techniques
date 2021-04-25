% PROGRAMMING EXERCISE 9: TIME-STEPPING WITH BDF
% Bruno Degli Esposti, Xingyu Xu
% 07/01/20 - 14/01/20
% Code tested in MATLAB only

%% Parameters

% of the PDE
c1 = 1;
c2 = 1;
coeff_a = @(x,y) c1;
coeff_r = @(x,y) c2;
uD = @(x,y,t) 0;
u0 = @(x,y) sin(pi*x) .* sin(pi*y);
u_ex = @(x,y,t) exp(-t) .* sin(pi*x) .* sin(pi*y);
f = @(x,y,t) (-1 + c1*2*pi^2 + c2) * u_ex(x,y,t);

% of the space discretization
dx = 0.1;
N = round(1/dx);
BC = [1;1;1;1];
[coord,elemNodeTable,boundary] = gen_mesh_rectangle(N,N,0,1,0,1,BC);
n_vertices = size(coord,1);

% of the time discretization
dt = 0.1;
T = 2;
nsteps = round(T/dt);
steporder = 1;

%% Solution

uh = BDF(coeff_a, coeff_r, f, uD, u0, ...
    coord, elemNodeTable, boundary, dt, T, steporder);

%% Validation of the results

u = zeros(n_vertices, nsteps+1);
for k = 0:nsteps
    u(:,k+1) = u_ex(coord(:,1),coord(:,2),dt*k);
end

e = uh-u;
for k = 0:nsteps
    t = dt*k;
    fprintf('At t = %f, err max = %f\n', t, norm(e(:,k+1),Inf));
end

% for k = 0:nsteps
%     trisurf(elemNodeTable,coord(:,1),coord(:,2),e(:,k+1));
%     shading interp; drawnow; pause;
% end



