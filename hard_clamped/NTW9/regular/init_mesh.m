function [vertices2, mesh2] = init_mesh(n)
gm = [3; 4; 0; 1; 1; 0; 0; 0; 1; 1];
model = createpde(1);
geometryFromEdges(model, decsg(gm));
mesh = generateMesh(model, 'Hmax', 1/n, 'GeometricOrder', 'linear');
Nv = size(mesh.Nodes, 2);
mesh = generateMesh(model, 'Hmax', 1/n);
[p, ~, t] = meshToPet(mesh);
Ne = size(p, 2) - Nv;
Nt = size(t, 2);
vertices = zeros(Nv+2*Ne, 2);
mesh = zeros(Nt, 9);
vertices(1: Nv+Ne, :) = p';
mesh(:, 1: 6) = t([1: 3, 5, 6, 4], :)';
vertices(Nv+Ne+1: Nv+2*Ne, :) = p(:, Nv+1: Nv+Ne)';
mesh(:, 7: 9) = t([5, 6, 4], :)' + Ne;
N = 2 * size(vertices, 1);
vertices2 = zeros(N, 2);
vertices2(1: 2: N, :) = vertices;
vertices2(2: 2: N, :) = vertices;
mesh2 = zeros(Nt, 18);
mesh2(:, 1: 2: 18) = 2 * mesh - 1;
mesh2(:, 2: 2: 18) = 2 * mesh;
end