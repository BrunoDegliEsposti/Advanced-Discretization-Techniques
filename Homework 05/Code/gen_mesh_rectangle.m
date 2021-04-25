% Expected input:
% N - number of subrectangles in x-direction
% M - number of subrectangles in y-direction 
% xmin, xmax, ymin, ymax - intervals containing the rectangle 
% sides = [S; E; N; W] - 4x1 vector specifying the boundary condition at the 
%                        South, East, North & West sides of the domain 
%                        (i.e. 1 = Dirichlet, 2 = Neumann, ...)

function [coord,elemNodeTable,boundary] = gen_mesh_rectangle(N, M, xmin, xmax, ymin, ymax, sides, el_type)
    if (el_type == 1)
        %Allocating memory for coord and elemNodeTable
        coord = zeros((N+1)*(M+1),2);
        elemNodeTable = zeros(2*M*N,3);
        
        % first the vertex coordinates
        deltax = (xmax - xmin)/N;
        deltay = (ymax - ymin)/M;
        for i = 0 : M
            for j = 0 : N
                coord(i*(N+1)+(j+1),:) = [xmin + j*deltax, ymin + i*deltay];
            end
        end
        
        % now the description of the elements
        count = 1;
        for i = 0 : M-1
            for j = 1 : N
                elemNodeTable(count,:) = [i*(N+1)+j , (i+1)*(N+1)+j+1, (i+1)*(N+1)+j];
                elemNodeTable(count+1,:) = [(i+1)*(N+1)+j+1, i*(N+1)+j , i*(N+1)+j+1];
                count = count+2;
            end
        end
        
        % finally the boundary segments with their respective flags
        boundary = zeros(2*N+2*M,3);
        count = 1;
        
        % at the base (y = ymin)
        for j = 1:N
            boundary(count,:) = [sides(1,1), j,j+1];
            count = count+1;
        end
        % at the lateral sides (x = xmin  or  x = xmax)
        for i = 1 : M
            boundary(count,:) = [sides(2,1),i*(N+1),(i+1)*(N+1)];
            boundary(count+1,:) = [sides(4,1),(i-1)*(N+1)+1,i*(N+1)+1];
            count = count+2;
        end
        % at the top (y = ymax)
        for j = 1:N
            boundary(count,:) = [sides(3,1),M*(N+1)+j,M*(N+1)+j+1];
            count = count+1;
        end
    elseif (el_type == 2)
        coord = zeros((2*N+1)*(2*M+1),2);
        elemNodeTable = zeros(2*M*N,6);
        boundary = zeros(2*N+2*M,4);
        deltax = (xmax - xmin)/N;
        deltay = (ymax - ymin)/M;
        
        % first the vertex coordinates
        for i = 0:M
            for j = 0:N
                coord(i*(N+1)+j+1,:) = [xmin + j*deltax, ymin + i*deltay];
            end
        end
        ecn = (N+1)*(M+1); % end of corner nodes
        for i = 0:M
            for j = 0:N-1
                coord(ecn+i*N+j+1,:) = [xmin + j*deltax + deltax/2, ymin + i*deltay];
            end
        end
        ehem = ecn + N*(M+1); % end of horizontal edges midpoints
        for i = 0:M-1
            for j = 0:2*N
                coord(ehem+i*(2*N+1)+j+1,:) = [xmin + j*(deltax/2), ymin + i*deltay + deltay/2];
            end
        end
        
        % now the description of the elements
        count = 1;
        for i = 0:M-1
            for j = 0:N-1
                U = [1,N+3,N+2, ecn+N+1,ehem+1,ehem+2]; % first triangle, bottom left
                hincU = [1,1,1,1,2,2]; % horizontal increment for the indices
                vincU = [N+1,N+1,N+1, N,2*N+1,2*N+1]; % vertical increment for the indices
                
                L = [1,2,N+3, ehem+3,ehem+2,ecn+1]; % second triangle, bottom left
                hincL = [1,1,1,2,2,1]; % horizontal increment for the indices
                vincL = [N+1,N+1,N+1, 2*N+1,2*N+1,N]; % vertical increment for the indices
                
                elemNodeTable(count,:) = U + i*vincU + j*hincU;
                elemNodeTable(count+1,:) = L + i*vincL + j*hincL;
                count = count + 2;
            end
        end
        
        % finally the boundary segments with their respective flags
        count = 1;
        
        % at the base (y = ymin)
        for j = 1:N
            boundary(count,:) = [sides(1,1), j, ecn+j, j+1];
            count = count+1;
        end
        % at the lateral sides (x = xmin  or  x = xmax)
        for i = 1:M
            boundary(count,:) = [sides(2,1), i*(N+1), ehem+i*(2*N+1), (i+1)*(N+1)];
            boundary(count+1,:) = [sides(4,1), (i-1)*(N+1)+1, ehem+(i-1)*(2*N+1)+1, i*(N+1)+1];
            count = count+2;
        end
        % at the top (y = ymax)
        for j = 1:N
            boundary(count,:) = [sides(3,1), M*(N+1)+j, ecn+M*N+j, M*(N+1)+j+1];
            count = count+1;
        end
    else
        error('Invalid element type');
    end
end















