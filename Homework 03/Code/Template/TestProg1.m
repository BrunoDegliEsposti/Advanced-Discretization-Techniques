%Tests AssembleSystem.m in the case of a stationary elliptic problem of the format

%  - nabla \cdot (a nabla u) + ru  = f      in Omega
%                                u = u_D    on Gamma

% Luca Wester <luca.wester@fau.de>, July 2019 

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
f = @(x,y) -(d1a(x,y).*grad_u_ex_1(x,y) + d2a(x,y).*grad_u_ex_2(x,y) + coeff_a(x,y).*lap_u(x,y)) + coeff_r(x,y).*u_ex(x,y);

%Generate mesh on [0,1] x [0,2]
[coord,elemNodeTable,boundary] = gen_mesh_rectangle(40,40,0,1,0,2,[1;1;1;1]);

%Get Dirichlet & free node indices
DirichletNodes = unique(boundary(find(boundary(:,1)==1),2:3));
FreeNodes = setdiff(1:size(coord,1),DirichletNodes);

%Initialize uh with Dirichlet values
uh = zeros(size(coord,1),1);
uh(DirichletNodes,1) = u_ex(coord(DirichletNodes,1),coord(DirichletNodes,2));

%Generate system matrices & right hand side
[A,M] = AssembleMatrices(coord,elemNodeTable,coeff_a,coeff_r);
[rhs] = AssembleRHS(coord,elemNodeTable,f);

%Update rhs with Dirichlet values
rhs = rhs - (A+M)*uh;

%Solve on only the free nodes
uh(FreeNodes,1) = (A(FreeNodes,FreeNodes)+M(FreeNodes,FreeNodes))\rhs(FreeNodes,1);

%% at this point 'uh' contains the solution at each vertex
%% we plot it with
u = u_ex(coord(:,1),coord(:,2));
subplot(1,2,1), trisurf(elemNodeTable, coord(:,1), coord(:,2), uh), shading interp;
subplot(1,2,2), trisurf(elemNodeTable, coord(:,1), coord(:,2), u), shading interp;
%
figure(2)
trisurf(elemNodeTable,coord(:,1),coord(:,2),u-uh), shading interp;













