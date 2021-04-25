coeff_a = @(x,y) x+sin(y).^2;
coeff_r = @(x,y) 2.*x+2.*y;
N = 3;
[coord,elemNodeTable,boundary] = gen_mesh_rectangle(N, N, ...
        0, 1, 0, 2, [1;1;1;1], 2);
[A,M] = AssembleMatrices(coord,elemNodeTable,coeff_a,coeff_r);
[Aref,Mref] = AssembleMatricesRef(coord,elemNodeTable,coeff_a,coeff_r);
normest(A-Aref,'inf')