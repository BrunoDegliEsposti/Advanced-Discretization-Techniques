%Tests AssembleSystem.m in the case of a stationary elliptic problem of the format

%  - nabla \cdot (a nabla u) + ru  = f      in Omega
%                                u = u_D    on Gamma

% ############################## Test case #################################### 
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

% #############################################################################


%TODO: Load mesh, assemble system, introduce Dirichlet b.c., solve


% #############################################################################


%% at this point 'uh' contains the solution at each vertex
%% we plot it with
u = u_ex(coord(:,1),coord(:,2));
subplot(1,2,1), trisurf(elemNodeTable, coord(:,1), coord(:,2), uh), shading interp;
subplot(1,2,2), trisurf(elemNodeTable, coord(:,1), coord(:,2), u), shading interp;
%
figure(2)
trisurf(elemNodeTable,coord(:,1),coord(:,2),u-uh), shading interp;
