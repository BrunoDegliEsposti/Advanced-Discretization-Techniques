classdef RTMesh < handle
    % RTMesh Container holding all information about the mesh
    %   Reads *.msh files and builds all the datastructures needed for the 
    %   Raviart-Thomas finite element.
    %   To save storage space and to identify objects uniquely, relations
    %   are stored using indices, e.g. for each triangle we store the
    %   indices of its three vertices, the coordinates themselves are
    %   stored in another list, see the property descriptions below.
    %   Utility functions are provided to query e.g. geometric information
    %   about edges or triangles.
	properties 
        % coordinates - List of (x,y) coordinates of mesh points 
		coordinates
        
        % triangles - List of triangles
        %   Each triangle consists of the indices (!) of the three
        %   coordinates of its vertices. 
        %   Example: To access the coordinates of vertex 2 of triangle 15:
        %       mesh.coordinates(mesh.triangles(15, 2), :)
		triangles
        
        % nodes2edge - Map between nodes and edges
        %   This is a sparse matrix with nodes2edge(i, j) = k if the nodes
        %   i and j belong to edge k (k = 0 if (i, j) is no edge).
		nodes2edge
        
        % edge2triangle - Matrix containing nodes and adjacent triangles
        %   A (numEdges x 4)-matrix, where the first two entries in each
        %   row are the vertices of the edge, the second two entries the
        %   indices of adjacent triangles. The first triangle corresponds
        %   to K+ (i.e. here the orientation is encoded), the second
        %   triangle to K- and may be 0 (for boundary edges).
		edge2triangle
        
        % numEdges - Number of edges in the mesh
		numEdges
        
        % numTriangles - Number of triangles in the mesh
        numTriangles
        
        % boundaryEdges - List of edges on the boundary
        %   The format is the same as in nodes2edge, except that we don't
        %   store the second triangle the edge belongs to (for obious 
        %   reasons).
		boundaryEdges
	end

	methods
		function mesh = RTMesh(filepath)
            loadMsh(mesh, filepath);
            initConnectivities(mesh);
        end
        
        function point = calcMidpoint(mesh, vertices)
            point = mean(mesh.coordinates(vertices, :));
        end
        
        function area = calcTriangleArea(mesh, triangle)
            % calcTriangleArea - calculates area of given triangle
            %   triangle: List of three vertex indices of triangle
            area = 0.5 * det([1, 1, 1 ; mesh.coordinates(triangle, :)']);
        end
        
        function length = calcEdgeLength(mesh, edge)
            % calcEdgeLength - calculates length of given edge
            %   edge: The two indices of the edge
            length = norm(mesh.coordinates(edge(1),:) - mesh.coordinates(edge(2),:));
        end
        
        function index = findEdgeIndex(mesh, nodes)
            index = mesh.nodes2edge(nodes(1), nodes(2));
        end
        
        function plotConstantFunction(mesh, coefficients, plotTitle)
            % plotConstantFunction  Plots piecewise constant function
            %   coefficients: Array of size numFaces, each entry being the
            %   coefficient of the Lagrange0 basis function corresponding
            %   to the triangle
            %   plotTitle: Title of the plot window
            
            trias = zeros(mesh.numTriangles, 3);
            coords = zeros(mesh.numTriangles * 3, 2);
            values = repmat(coefficients', 3, 1);

            for i=1:mesh.numTriangles
                trias(i,:) = 3 * (i - 1) + (1:3);
                coords(trias(i,:), :) = mesh.coordinates(mesh.triangles(i,:), :);
            end

            trisurf(trias, coords(:, 1), coords(:, 2), values(:));
            shading interp;
            title(plotTitle)
            view(-60,50);
        end
        
        function plotRTFunction(mesh, coefficients, plotTitle)
            % plotRTFunction  Plots piecewise constant function
            %   coefficients: Array of size numEdges, each entry being the
            %   coefficient of the RT basis function associated with that
            %   edge
            %   plotTitle: Title of the plot window
            
            trias = zeros(mesh.numTriangles, 3);
            coords = zeros(mesh.numTriangles * 3, 2);
            values = zeros(mesh.numTriangles * 3, 2);

            for i=1:mesh.numTriangles
                trias(i,:) = 3 * (i - 1) + (1:3);
                coords(trias(i,:), :) = mesh.coordinates(mesh.triangles(i,:), :);
                
                faces = mesh.nodes2edge(mesh.triangles(i, [2 3 1]), mesh.triangles(i, [3 1 2]));
                faceNumbers = diag(faces);
                
                area = mesh.calcTriangleArea(mesh.triangles(i, :));

                % evaluate the Raviart-Thomas basis functions
                for k = 1:3
                    vertex = mesh.coordinates(mesh.triangles(i, k), :);
                    
                    for l = 1:3
                        edgeOppositeVertex = mesh.coordinates(mesh.triangles(i, l), :);
                        signum = 1 - 2 * (i == mesh.edge2triangle(faceNumbers(l), 4));
                        
                        edgeLength = mesh.calcEdgeLength(mesh.edge2triangle(faceNumbers(l), 1:2));
                        values(3 * (i-1) + k, :) = values(3 * (i-1) + k, :) + signum * edgeLength / (2 * area) * (vertex - edgeOppositeVertex) * coefficients(faceNumbers(l));
                    end
                end
            end

            subplot(2,1,1)
            trisurf(trias, coords(:, 1), coords(:, 2), values(:, 1));
            shading interp;
            title([plotTitle , ' (x)']);
            view(-60,50);
            
            subplot(2,1,2)
            trisurf(trias, coords(:, 1), coords(:, 2), values(:, 2));
            shading interp;
            title([plotTitle , ' (y)']);
            view(-60,50);
        end
    end
    
    % These are internal methods
    methods (Access=private)
        function loadMsh(mesh, filepath)
            % loadMsh  Reads .msh file
            %   We ignore the formal specification and just extract the 
            %   relevant information using hardcoded positions. This works
            %   as long as the file only contains the absolute minimum 
            %   information (e.g. no physical groups)
            
            numMshNodes = dlmread(filepath, '', [4 0 4 0]);
            numMshElements = dlmread(filepath, '', [7 + numMshNodes 0 7 + numMshNodes 0]);

            % extract coordinatess of vertices
            mesh.coordinates = dlmread(filepath, '', [5 1 4 + numMshNodes 2]);
            
            % mesh elements may be points, lines or triangles, we are only
            % interested in the latter (identified by second column == 2)
            mshElements = dlmread(filepath, '', [8 + numMshNodes 0 7 + numMshNodes + numMshElements 7]);

            mesh.triangles = mshElements((mshElements(:, 2) == 2), 6:8);
        end
        
        function initConnectivities(mesh)
            % initConnectivities  Build relation tables between mesh objects
            
            nodes2triangles=sparse(size(mesh.coordinates,1),size(mesh.coordinates,1));
            for j=1:size(mesh.triangles,1)
                nodes2triangles(mesh.triangles(j,:),mesh.triangles(j,[2 3 1]))=...
                  nodes2triangles(mesh.triangles(j,:),mesh.triangles(j,[2 3 1]))+j*eye(3,3);
            end

            B=nodes2triangles+nodes2triangles';
            [i,j]=find(triu(B));
            mesh.nodes2edge=sparse(i,j,1:size(i,1),size(mesh.coordinates,1),size(mesh.coordinates,1));
            mesh.nodes2edge=mesh.nodes2edge+mesh.nodes2edge';
            
            mesh.numEdges = size(i,1);
            mesh.numTriangles = size(mesh.triangles, 1);
            
            % to generate element of edge
            mesh.edge2triangle=zeros(size(i,1),4);
            for m = 1:size(mesh.triangles,1)
                for k = 1:3
                    p = mesh.nodes2edge(mesh.triangles(m,k),mesh.triangles(m,rem(k,3)+1));
                    if mesh.edge2triangle(p,1)==0
                        mesh.edge2triangle(p,:)=[mesh.triangles(m,k) mesh.triangles(m,rem(k,3)+1)  nodes2triangles(mesh.triangles(m,k),mesh.triangles(m,rem(k,3)+1)) ...
                            nodes2triangles(mesh.triangles(m,rem(k,3)+1),mesh.triangles(m,k))];
                        
                    end
                end
            end

            % boundary edges are identified by their second triangle being 0
            mesh.boundaryEdges = mesh.edge2triangle(mesh.edge2triangle(:,4) == 0, 1:3);
        end
    end
end
