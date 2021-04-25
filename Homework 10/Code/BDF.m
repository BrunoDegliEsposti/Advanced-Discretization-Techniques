% PROGRAMMING EXERCISE 9: TIME-STEPPING WITH BDF
% Bruno Degli Esposti, Xingyu Xu
% 07/01/20 - 14/01/20
% Code tested in MATLAB only

function [uh] = BDF(coeff_a, coeff_r, f, uD, u0, ...
    coord, elemNodeTable, boundary, dt, T, steporder)
    %BDF Time-stepping with BDF for the Advection-Diffusion-Reaction PDE
    if (all(steporder ~= [1,2,3]))
        error('Invalid step order')
    end
    
    t = 0;
    n = 1; % defined so that U^n in the BDF formulas is uh(:,n)
    nsteps = round(T/dt);
    n_vertices = size(coord,1);
    boundary_vertices = unique([boundary(:,2);boundary(:,3)])';
    
    % initialize uh with the initial condition u0
    uh = zeros(n_vertices, nsteps+1);
    uh(:,1) = u0(coord(:,1),coord(:,2));
    
    % precompute stiffness matrix, reaction matrix, mass matrix
    [A,R] = AssembleMatrices(coord,elemNodeTable,coeff_a,coeff_r);
    [~,M] = AssembleMatrices(coord,elemNodeTable,@(x,y)0,@(x,y)1);
    
    % precompute Crank-Nicolson Matrices
    CNM_lhs = M + (dt/2)*A + (dt/2)*R;
    CNM_rhs = M - (dt/2)*A - (dt/2)*R;
    
    % precompute the BDF Matrix
    if (steporder == 1)
        BDFM = M + dt*A + dt*R;
    elseif (steporder == 2)
        BDFM = M + (2*dt/3)*A + (2*dt/3)*R;
    else
        BDFM = M + (6*dt/11)*A + (6*dt/11)*R;
    end
    
    % precompute the effects of Dirichlet BC on the matrices
    for v_idx = boundary_vertices
        CNM_lhs(v_idx,:) = 0;
        CNM_lhs(v_idx,v_idx) = 1;
        BDFM(v_idx,:) = 0;
        BDFM(v_idx,v_idx) = 1;
    end
    
    % precompute the time-dependent forcing term
    forcing = zeros(n_vertices, nsteps+1);
    for k = 0:nsteps
        forcing(:,k+1) = AssembleRHSt(coord,elemNodeTable,f,dt*k);
    end
    
    % preliminary Crank-Nicolson steps
    while (n < steporder)
        rhs = CNM_rhs*uh(:,n) + (dt/2)*(forcing(:,n+1)+forcing(:,n));
        for v_idx = boundary_vertices
            v = coord(v_idx,:);
            rhs(v_idx) = uD(v(1),v(2),t+dt);
        end
        uh(:,n+1) = CNM_lhs\rhs;
        t = t + dt;
        n = n + 1;
    end
    
    % BDF steps
    while (n <= nsteps)
        if (steporder == 1)
            rhs = M*uh(:,n) + dt*forcing(:,n+1);
        elseif (steporder == 2)
            rhs = (4/3)*M*uh(:,n) - (1/3)*M*uh(:,n-1) + ...
                (2/3)*dt*forcing(:,n+1);
        else
            rhs = (18/11)*M*uh(:,n) - (9/11)*M*uh(:,n-1) + ...
                (2/11)*M*uh(:,n-2) + (6/11)*dt*forcing(:,n+1);
        end
        % Dirichlet BC
        for v_idx = boundary_vertices
            v = coord(v_idx,:);
            rhs(v_idx) = uD(v(1),v(2),t+dt);
        end
        uh(:,n+1) = BDFM\rhs;
        t = t + dt;
        n = n + 1;
    end
end








