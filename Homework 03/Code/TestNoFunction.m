%% Parameters
c1 = 1;
c2 = 1;
coeff_a = @(x,y) c1;
coeff_r = @(x,y) c2;
u_ex = @(x,y,t) exp(-t) .* sin(pi*x) .* sin(pi*y);
f = @(x,y,t) (-1 + c1*2*pi^2 + c2) * u_ex(x,y,t);

N = 32;
deltax = 1/(N-1);
deltay = 1/(N-1);

T = 2;
deltat = 0.1;
nsteps = round(T/deltat);
theta = 0;

%% Assembly

BC = [1;1;1;1];
[coord,elemNodeTable,boundary] = gen_mesh_rectangle(N,N,0,1,0,1,BC);
n_vertices = size(coord,1);
n_boundary_edges = size(boundary,1);

[A,wM] = AssembleMatrices(coord,elemNodeTable,coeff_a,coeff_r);
M = unweightedMassMatrix(coord,elemNodeTable);
A = A + wM;
S = M + deltat*theta*A;
W = M - deltat*(1-theta)*A;
% With this notation, theta-method iterations will be:
% S u_{k+1} = W u_k + deltat (theta rhs_{k+1} + (1-theta) rhs_k)

for i = 1:n_boundary_edges
    if boundary(i,1) == 1
        vertex_id = boundary(i,2);
        S(vertex_id,:) = 0;
        S(vertex_id,vertex_id) = 1;
        vertex_id = boundary(i,3);
        S(vertex_id,:) = 0;
        S(vertex_id,vertex_id) = 1;
    end
end

%% Time stepping scheme

t = 0;
uh = zeros(n_vertices, nsteps+1);
rhs = zeros(n_vertices, nsteps+1);

uh(:,1) = u_ex(coord(:,1),coord(:,2),0);
for k = 1:nsteps+1
    rhs(:,k) = AssembleRHS(coord,elemNodeTable,f,deltat*(k-1));
end

for k = 1:nsteps
    r = W*uh(:,k) + deltat * (theta*rhs(:,k+1) + (1-theta)*rhs(:,k));
    for i = 1:n_boundary_edges
        if boundary(i,1) == 1
            vertex_id = boundary(i,2);
            vertex = coord(vertex_id,:);
            r(vertex_id) = u_ex(vertex(1),vertex(2),t+deltat);
            vertex_id = boundary(i,3);
            vertex = coord(vertex_id,:);
            r(vertex_id) = u_ex(vertex(1),vertex(2),t+deltat);
        end
    end
    uh(:,k+1) = S\r;
    t = t + deltat;
end

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



