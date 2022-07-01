function [vertices2, mesh2] = init_mesh(n)
gm = [3; 4; 0; 1; 1; 0; 0; 0; 1; 1];
model = createpde(1);
geometryFromEdges(model, decsg(gm));
mesh = generateMesh(model, 'Hmax', 1/n, 'GeometricOrder', 'linear');
Nv = size(mesh.Nodes, 2);
mesh = generateMesh(model, 'Hmax', 1/n);
[p, ~, t] = meshToPet(mesh);
%pdemesh(p,e,t, 'NodeLabels', 'on');
Ne = size(p, 2) - Nv;
Nt = size(t, 2);
vertices = zeros(Nv+2*Ne+Nt, 2);
mesh = zeros(Nt, 10);
vertices(1: Nv+Ne, :) = p';
mesh(:, 1: 6) = t([1: 3, 5, 6, 4], :)';
vertices(Nv+Ne+1: Nv+2*Ne, :) = p(:, Nv+1: Nv+Ne)';
mesh(:, 7: 9) = t([5, 6, 4], :)' + Ne;
vertices(Nv+2*Ne+1: Nv+2*Ne+Nt, :) = (p(:, t(1, :)) + p(:, t(2, :)) + p(:, t(3, :)))' / 3;
mesh(:, 10) = Nv+2*Ne+1: Nv+2*Ne+Nt;
N = 2 * size(vertices, 1);
vertices2 = zeros(N+Nv, 2);
vertices2(1: 2: N, :) = vertices;
vertices2(2: 2: N, :) = vertices;
mesh2 = zeros(Nt, 23);
mesh2(:, 1: 2: 20) = 2 * mesh - 1;
mesh2(:, 2: 2: 20) = 2 * mesh;
vertices2(N+1: N+Nv, :) = p(:, 1: Nv)';
mesh2(:, 21: 23) = t(1: 3, :)' + N;
end