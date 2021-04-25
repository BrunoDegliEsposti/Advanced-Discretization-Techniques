% PROGRAMMING EXERCISE 10: A FRACTIONAL THETA-SCHEME
% Bruno Degli Esposti, Xingyu Xu
% 14/01/20 - 21/01/20
% Code tested in MATLAB only

function [uh] = FractionalTheta(coeff_a, coeff_r, f, uD, u0, ...
    coord, elemNodeTable, boundary, dt, T, theta, alpha)
    %FRACTIONALTHETA Fractional theta stepping method for parabolic PDEs
    if (theta <= 0 || theta >= 1/2)
        error('Error: theta is not in (0,1/2)')
    end
    if (alpha <= 0 || alpha >= 1)
        error('Error: alpha is not in (0,1)')
    end
    theta_prime = 1-2*theta;
    
    t = 0;
    n = 1; % defined so that U^n in the formulas is stored in uh(:,n)
    nsteps = round(T/dt);
    n_vertices = size(coord,1);
    boundary_vertices = unique([boundary(:,2);boundary(:,3)])';
    
    % initialize uh with the initial condition u0
    uh = zeros(n_vertices, nsteps+1);
    uh(:,1) = u0(coord(:,1),coord(:,2));
    
    % precompute stiffness matrix, reaction matrix, mass matrix
    [A,R] = AssembleMatrices(coord,elemNodeTable,coeff_a,coeff_r);
    [~,M] = AssembleMatrices(coord,elemNodeTable,@(x,y)0,@(x,y)1);
    
    % precompute the Fractional Theta Method matrices
    FTM_lhs_1 = M + alpha*theta*dt*(A+R);
    FTM_rhs_1 = M - (1-alpha)*theta*dt*(A+R);
    FTM_lhs_2 = M + (1-alpha)*theta_prime*dt*(A+R);
    FTM_rhs_2 = M - alpha*theta_prime*dt*(A+R);
    FTM_lhs_3 = FTM_lhs_1;
    FTM_rhs_3 = FTM_rhs_1;
    
    % precompute the effects of Dirichlet BC on the matrices
    for v_idx = boundary_vertices
        FTM_lhs_1(v_idx,:) = 0;
        FTM_lhs_1(v_idx,v_idx) = 1;
        FTM_lhs_2(v_idx,:) = 0;
        FTM_lhs_2(v_idx,v_idx) = 1;
        FTM_lhs_3(v_idx,:) = 0;
        FTM_lhs_3(v_idx,v_idx) = 1;
    end
    
    % fractional theta stepping method.
    % naming convention: _0 means ^n, _1 means ^{n+\theta},
    % _2 means ^{n+1-\theta}, _3 means ^{n+1}
    while (n <= nsteps)
        forcing_0 = AssembleRHSt(coord,elemNodeTable,f,dt*(n-1));
        forcing_2 = AssembleRHSt(coord,elemNodeTable,f,dt*(n-theta));
        % First equation
        rhs_1 = FTM_rhs_1*uh(:,n) + theta*dt*forcing_0;
        for v_idx = boundary_vertices
            v = coord(v_idx,:);
            rhs_1(v_idx) = uD(v(1),v(2),t+dt*theta);
        end
        uh_1 = FTM_lhs_1 \ rhs_1;
        % Second equation
        rhs_2 = FTM_rhs_2*uh_1 + theta_prime*dt*forcing_2;
        for v_idx = boundary_vertices
            v = coord(v_idx,:);
            rhs_2(v_idx) = uD(v(1),v(2),t+dt*(1-theta));
        end
        uh_2 = FTM_lhs_2 \ rhs_2;
        % Third equation
        rhs_3 = FTM_rhs_3*uh_2 + theta*dt*forcing_2;
        for v_idx = boundary_vertices
            v = coord(v_idx,:);
            rhs_3(v_idx) = uD(v(1),v(2),t+dt);
        end
        uh(:,n+1) = FTM_lhs_3 \ rhs_3;
        t = t + dt;
        n = n + 1;
    end
end








