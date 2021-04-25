% PROGRAMMING EXERCISE 11: CONVECTION-DOMINATED PROBLEMS WITH SUPG
% Bruno Degli Esposti, Xingyu Xu
% 21/01/20 - 04/02/20
% Code tested in MATLAB only

%% Parameters

% of the PDE
a = 1;
r = 1;
coeff_a = @(x,y) a;
coeff_c = @(x,y) [3;-7];
coeff_r = @(x,y) r;
uD = @(x,y,t) 0;
u0 = @(x,y) sin(pi*x) .* sin(pi*y);
u_ex = @(x,y,t) exp(-t) .* sin(pi*x) .* sin(pi*y);
f = @(x,y,t) (-1 + 2*pi^2*a + r) * u_ex(x,y,t) + ...
    pi*exp(-t)*dot(coeff_c(x,y), [cos(pi*x).*sin(pi*y); sin(pi*x).*cos(pi*y)]);

% of the space discretization
dx = 0.05;
N = round(1/dx);
BC = [1;1;1;1];
[coord,elemNodeTable,boundary] = gen_mesh_rectangle(N,N,0,1,0,1,BC);
n_vertices = size(coord,1);

% of the time discretization
dt = 0.1;
T = 2;
nsteps = round(T/dt);

%% Solution

uh = CrankNicolson(coeff_a, coeff_c, coeff_r, f, uD, u0, ...
    coord, elemNodeTable, boundary, dt, T);

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



