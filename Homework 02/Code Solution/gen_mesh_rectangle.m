% generates a uniform mesh of the rectangle
%        [xmin, xmax] x [ymin, ymax]
% dividing it into  NxM  subrectangles
% and then each subrectangle into two triangles

% the numbering of vertices and elements is as follows
%
% (N+1)M +---+---+---+---+(N+1)(M+1)
%        | / | / | / | / |
%        +---+---+---+---+
%        | / | / | / | / |
%    N+2 +---+---+---+---+N+2+N
%        |1/2|3/4| / | / |
%        +---+---+---+---+
%        1   2   ...     N+1

% Luca Wester <luca.wester@fau.de>, July 2019 
% based on finite element code for a 2d diffusion-reaction 
% system by Pedro Morin, Universidat Nacional del Litoral  
% Santa Fe, Argentinia 

% Expected input:
% N - number of subrectangles in x-direction
% M - number of subrectangles in y-direction 
% xmin, xmax, ymin, ymax - intervals containing the rectangle 
% sides = [S; E; N; W] - 4x1 vector specifying the boundary condition at the 
%                        South, East, North & West sides of the domain 
%                        (i.e. 1 = Dirichlet, 2 = Neumann, ...)

function [coord,elemNodeTable,boundary] = gen_mesh_rectangle(N, M, xmin, xmax, ymin, ymax,sides)


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


end 