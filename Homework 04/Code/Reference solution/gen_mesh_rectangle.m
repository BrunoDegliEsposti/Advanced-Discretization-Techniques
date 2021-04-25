% generates a uniform mesh of the rectangle
%        [xmin, xmax] x [ymin, ymax]
% dividing it into  NxM  subrectangles
% and then each subrectangle into two triangles

% the numbering of vertices and elements is as follows
%
% Linear Lagrange:
%
% (N+1)M +---+---+---+---+(N+1)(M+1)
%        | / | / | / | / |
%        +---+---+---+---+
%        | / | / | / | / |
%    N+2 +---+---+---+---+N+2+N
%        |1/2|3/4| / | / |
%        +---+---+---+---+
%        1   2   ...     N+1
%
%
% Quadratic Lagrange:
%
%        .                 .
%        .          .      .         .       .
%        .        .        .       .         .
%        |      .          |     .           .
%        |                 |                 | 
%        |      R+N+1      |      R+N+2      |                 .
%        +-------*---------+--------*--------+---- . . .       .
%        |              /  |              /  |                 .
%        |    1       /    |     3      /    |                 |
%        |         /       |         /       |       ....      |
%       S+1      S+2      S+3     S+4        |                 |
%        |     /       2   |     /      4    |                 |
%        |  /              |  /              |                 |
%        +-------*---------+-------*---------+-----------------+
%        1      R+1        2      R+2       ...               N+1
%        
%
% where R = (M+1)(N+1) = number of total triangle corner points (which are marked with +)
% and S = R + (M+1)N = R + (number of points marked with *)
%
% Numbering the corner points first ensures that the corresponding P1-mesh can 
% be extracted from the first 3 columns of elemNodeTable
%
% Luca Wester <luca.wester@fau.de>, November 2019 
% based on finite element code for a 2d diffusion-reaction     %
% system by Pedro Morin, Universidat Nacional del Litoral       %
% Santa Fe, Argentinia 

% Expected input:
% N - number of subrectangles in x-direction
% M - number of subrectangles in y-direction 
% xmin, xmax, ymin, ymax - intervals containing the rectangle 
% sides = [S; E; N; W] - 4x1 vector specifying the boundary condition at the 
%                        South, East, North & West sides of the domain 
%                        (i.e. 1 = Dirichlet, 2 = Neumann, ...)
% el_type - element type (i.e. 1 = linear, 2 = quadratic, ...)

function [coord,elemNodeTable,boundary] = gen_mesh_rectangle(N, M, xmin, xmax, ymin, ymax,sides,el_type)


%Allocating memory for coord and elemNodeTable
if (el_type == 1)
  
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
  
  % first the corner nodes
  deltax = (xmax - xmin)/N;
  deltay = (ymax - ymin)/M;
  for i = 0 : (M)
    for j = 0 : (N)
      coord(i*(N+1)+(j+1),:) = [xmin + j*deltax, ymin + i*deltay];
    end
  end
  
  % then the horizontal midpoints of the corner nodes
  R = (M+1)*(N+1);
  for i = 0:M
    for j = 1:N
      coord((M+1)*(N+1)+i*N+j,:) = [xmin + deltax/2 + (j-1)*deltax,ymin + i*deltay];
    end
  end
  
  % finally the rows lying vertically between the already listed points
  S = R + (M+1)*N;
  for i = 0:(M-1)
    for j = 0:(2*N)
      coord((M+1)*(N+1)+(M+1)*N+i*(2*N+1)+j+1,:) = [xmin + j/2*deltax,ymin+deltay/2+i*deltay];
    end
  end


  % now the description of the elements
  count = 1;
  for i = 1 : M
    for j = 1 : N
      elemNodeTable(count,:) = [(i-1)*(N+1)+j, i*(N+1)+j+1, i*(N+1)+j, ...
                                R+i*N+j, S+(i-1)*(2*N+1)+2*j-1, ...
                                S+(i-1)*(2*N+1)+2*j];
      elemNodeTable(count+1,:) = [(i-1)*(N+1)+j,(i-1)*(N+1)+j+1,i*(N+1)+j+1, ...
                                  S+(i-1)*(2*N+1)+2*j+1,...
                                  S+(i-1)*(2*N+1)+2*j,R+(i-1)*N+j];
      count = count+2;
    end
  end
  
  % finally the boundary segments with their respective flags  
  boundary = zeros(2*N+2*M,4);
  count = 1;
  
  % at the base (y = ymin)
  for j = 1:N
    boundary(count,:) = [sides(1,1),j,R+j,j+1];
    count = count+1;
  end
  % at the lateral sides (x = xmin  or  x = xmax)
  for i = 1 : M
      boundary(count,:) = [sides(2,1),i*(N+1),S+i*(2*N+1),(i+1)*(N+1)];
      boundary(count+1,:) = [sides(4,1),1+(i-1)*(N+1),S+1+(i-1)*(2*N+1),1+i*(N+1)];
      count = count+2;
  end
  % at the top (y = ymax)
  for j = 1:N
    boundary(count,:) = [sides(3,1),M*(N+1)+j,R+(M*N)+j,M*(N+1)+j+1];
    count = count+1;
  end

else 

  error('Element type not supported!')

end


end 