% PROGRAMMING EXERCISE 6: UZAWA’S ALGORITHM
% Bruno Degli Esposti, Xingyu Xu
% 3/12/19 - 10/12/19
% Code tested in MATLAB only

function [u,p,k] = uzawa(A,B,M,f,g,p0,omega,tol,kmax)
    u = zeros(size(A,1),1);
    p = p0;
    pold = p0;
    for k = 1:kmax
        u = A\(f-B'*p);
        t = -omega*(g-B*u);
        p = pold + M\t;
        if dot(p-pold,t) < tol*tol
            break;
        end
        pold = p;
    end
end