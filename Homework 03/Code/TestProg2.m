%% Parameters

% of the differential equation
c1 = 1;
c2 = 1;
coeff_a = @(x,y) c1;
coeff_r = @(x,y) c2;
u_ex = @(x,y,t) exp(-t) .* sin(pi*x) .* sin(pi*y);
f = @(x,y,t) (-1 + c1*2*pi^2 + c2) * u_ex(x,y,t);

% of the space discretization
N = 32;
deltax = 1/(N-1);
deltay = 1/(N-1);
BC = [1;1;1;1];
[coord,elemNodeTable,boundary] = gen_mesh_rectangle(N,N,0,1,0,1,BC);
n_vertices = size(coord,1);

% of the time discretization
deltat = 0.1;
T = 2;
nsteps = round(T/deltat);
theta = 0.5;

%% Solution

uh = TimeStepping(coeff_a, coeff_r, u_ex, f, coord, elemNodeTable, boundary, deltat, T, theta);

%% Validation of the results

u = zeros(n_vertices, nsteps+1);
for k = 1:nsteps+1
    u(:,k) = u_ex(coord(:,1),coord(:,2),deltat*(k-1));
end

e = uh-u;
for k = 1:nsteps+1
    t = deltat*(k-1);
    fprintf('At t = %f, err max = %f\n', t, norm(e(:,k),Inf));
end

% for k = 1:nsteps+1
%     trisurf(elemNodeTable,coord(:,1),coord(:,2),e(:,k));
%     shading interp; drawnow; pause;
% end



