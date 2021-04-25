% PROGRAMMING EXERCISE 7: RAVIART-THOMAS ELEMENTS POINT A)
% Bruno Degli Esposti, Xingyu Xu
% 10/12/19 - 17/12/19
% Code tested in MATLAB only

% Data: exact solution u and right-hand side f
u = @(x,y) sin(pi*x).*sin(pi*y);
f = @(x,y) 2*pi^2*sin(pi*x).*sin(pi*y);

% Load mesh from file
mesh = RTMesh('meshes/square-3.msh');

numRTDoFs = mesh.numEdges;
numPDoFs = mesh.numTriangles;

% Implicitly we assume the following order of degrees of freedoms in our
% linear systems: first all RT dofs, then all P0 dofs. The RT dofs are
% numbered according to the edge number (i.e. the number in mesh.nodes2edge)
% the number of the P0 dofs matches the numbering of mesh.triangles.

% Assemble matrix blocks A and B corresponding to the bilinear forms a and b
A = sparse(numRTDoFs,numRTDoFs);
B = sparse(numPDoFs,numRTDoFs);
M = [2,0,1,0,1,0;
     0,2,0,1,0,1;
     1,0,2,0,1,0;
     0,1,0,2,0,1;
     1,0,1,0,2,0;
     0,1,0,1,0,2];

for j = 1:mesh.numTriangles
    v1 = mesh.triangles(j,1);
    v2 = mesh.triangles(j,2);
    v3 = mesh.triangles(j,3);
    x1 = mesh.coordinates(v1,:)';
    x2 = mesh.coordinates(v2,:)';
    x3 = mesh.coordinates(v3,:)';
    e1 = mesh.nodes2edge(v2,v3);
    e2 = mesh.nodes2edge(v3,v1);
    e3 = mesh.nodes2edge(v1,v2);
    sigma1 = (-1)^(j~=mesh.edge2triangle(e1,3));
    sigma2 = (-1)^(j~=mesh.edge2triangle(e2,3));
    sigma3 = (-1)^(j~=mesh.edge2triangle(e3,3));
    l1 = mesh.calcEdgeLength([v2,v3]);
    l2 = mesh.calcEdgeLength([v3,v1]);
    l3 = mesh.calcEdgeLength([v1,v2]);
    el_B = [sigma1*l1, sigma2*l2, sigma3*l3];
    B(j,[e1,e2,e3]) = B(j,[e1,e2,e3]) + el_B;
    c = 1/(48*abs(mesh.calcTriangleArea([v1,v2,v3])));
    N = [x1-x1,x1-x2,x1-x3; x2-x1,x2-x2,x2-x3; x3-x1,x3-x2,x3-x3];
    el_A = c*diag(el_B)*N'*M*N*diag(el_B);
    A([e1,e2,e3],[e1,e2,e3]) = A([e1,e2,e3],[e1,e2,e3]) + el_A;
end

% Assemble global stiffness matrix S
S = [A,B'; B,zeros(numPDoFs)];

% Assemble right-hand side terms stemming from forcing f
b = zeros(numRTDoFs+numPDoFs, 1);
for j = 1:mesh.numTriangles
    v = mesh.triangles(j,1:3);
    x = mesh.calcMidpoint(v);
    area = abs(mesh.calcTriangleArea(v));
    b(numRTDoFs+j) = -area*f(x(1),x(2));
end

% Incorporate Dirichlet conditions
% u_D = 0 on all the boundary, nothing to do

% Solve
solution = S\b;
sigmah = solution(1:numRTDoFs);
uh = solution(numRTDoFs+1:end);

% Plot results
figure(1)
mesh.plotConstantFunction(uh, 'u_h');

figure(2)
mesh.plotRTFunction(sigmah, 'Du_h');




