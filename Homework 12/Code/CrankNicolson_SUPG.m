% PROGRAMMING EXERCISE 11: CONVECTION-DOMINATED PROBLEMS WITH SUPG
% Bruno Degli Esposti, Xingyu Xu
% 21/01/20 - 04/02/20
% Code tested in MATLAB only

function [uh] = CrankNicolson_SUPG(coeff_a, coeff_c, coeff_r, f, ...
        uD, u0, coord, elemNodeTable, boundary, deltaT_SUPG, dt, T)
    %CRANKNICOLSONSUPG Crank-Nicolson scheme + Streamline Upwind
    % Petrov-Galerkin for the Advection-Diffusion-Reaction equation
    % with Dirichlet boundary conditions
    
    t = 0;
    n = 1; % defined so that U^n in the method is uh(:,n)
    nsteps = round(T/dt);
    n_vertices = size(coord,1);
    boundary_vertices = unique([boundary(:,2);boundary(:,3)])';
    
    % initialize uh with the initial condition u0
    uh = zeros(n_vertices, nsteps+1);
    uh(:,1) = u0(coord(:,1),coord(:,2));
    
    % precompute stiffness matrix, advection matrix,
    % reaction matrix, mass matrix, SUPG matrix
    [A,B,R] = AssembleMatrices(coord,elemNodeTable, ...
        coeff_a,coeff_c,coeff_r);
    [~,~,M] = AssembleMatrices(coord,elemNodeTable, ...
        @(x,y)0,@(x,y)[0;0],@(x,y)1);
    A_SUPG = AssembleMatrix_SUPG(coord,elemNodeTable, ...
        deltaT_SUPG,coeff_c,coeff_r);
    
    % precompute Crank-Nicolson Matrices
    CNM_lhs = M + (dt/2)*(A+B+R+A_SUPG);
    CNM_rhs = M - (dt/2)*(A+B+R+A_SUPG);
    
    % precompute the effects of Dirichlet BC on the matrices
    for v_idx = boundary_vertices
        CNM_lhs(v_idx,:) = 0;
        CNM_lhs(v_idx,v_idx) = 1;
    end
    
    % precompute the time-dependent forcing term + the SUPG term
    forcing = zeros(n_vertices, nsteps+1);
    for k = 0:nsteps
        forcing(:,k+1) = AssembleRHSt(coord,elemNodeTable,f,dt*k) ...
                       + AssembleRHSt_SUPG(coord,elemNodeTable, ...
                                           deltaT_SUPG,coeff_c,f,dt*k);
    end
    
    % Crank-Nicolson steps
    while (n <= nsteps)
        rhs = CNM_rhs*uh(:,n) + (dt/2)*(forcing(:,n+1)+forcing(:,n));
        % Dirichlet BC
        for v_idx = boundary_vertices
            v = coord(v_idx,:);
            rhs(v_idx) = uD(v(1),v(2),t+dt);
        end
        uh(:,n+1) = CNM_lhs\rhs;
        t = t + dt;
        n = n + 1;
    end
end








