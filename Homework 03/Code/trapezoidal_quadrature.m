function [integral] = trapezoidal_quadrature(vertices, area, f)
% vertices is a 3x2 matrix. Each row contains the coordinates
% of one vertex.
    v1 = vertices(1,:);
    v2 = vertices(2,:);
    v3 = vertices(3,:);
    m12 = 0.5*(v1+v2);
    m23 = 0.5*(v2+v3);
    m31 = 0.5*(v3+v1);
    integral = (area/3)*(f(m12(1),m12(2)) ...
                        +f(m23(1),m23(2)) ...
                        +f(m31(1),m31(2)));
end

