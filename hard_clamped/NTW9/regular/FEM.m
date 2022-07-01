function error = FEM(n, iota, lambda, mu, ux, uy, uxx, uxy, uyy, f, phi, phix, phiy, phixx, phixy, phiyy, T)
wg = [0.050844906370207; 0.050844906370207; 0.050844906370207; 0.116786275726379; 0.116786275726379; 0.116786275726379; 0.082851075618374; 0.082851075618374; 0.082851075618374; 0.082851075618374; 0.082851075618374; 0.082851075618374];
ld = [
    0.873821971016996, 0.063089014491502, 0.063089014491502;
    0.063089014491502, 0.873821971016996, 0.063089014491502;
    0.063089014491502, 0.063089014491502, 0.873821971016996;
    0.501426509658179, 0.249286745170910, 0.249286745170910;
    0.249286745170910, 0.501426509658179, 0.249286745170910;
    0.249286745170910, 0.249286745170910, 0.501426509658179;
    0.636502499121399, 0.310352451033785, 0.053145049844816;
    0.636502499121399, 0.053145049844816, 0.310352451033785;
    0.310352451033785, 0.636502499121399, 0.053145049844816
    0.310352451033785, 0.053145049844816, 0.636502499121399;
    0.053145049844816, 0.636502499121399, 0.310352451033785
    0.053145049844816, 0.310352451033785, 0.636502499121399;
];
[vertices, mesh] = init_mesh(n);
not_bdr = prod(vertices .* (1 - vertices), 2);
free = not_bdr ~= 0;
Nv = size(vertices, 1);
Nt = size(mesh, 1);
a = @(ux, uy, uxx, uxy, uyy, vx, vy, vxx, vxy, vyy) lambda * (ux(1) + uy(2)) * (vx(1) + vy(2)) + 2 * mu * (ux(1) * vx(1) + uy(2) * vy(2)) + mu * (uy(1) + ux(2)) * (vy(1) + vx(2)) + lambda * iota^2 * ((uxx(1) + uxy(2)) * (vxx(1) + vxy(2)) + (uxy(1) + uyy(2)) * (vxy(1) + vyy(2))) + 2 * mu * iota^2 * (uxx(1) * vxx(1) + uxy(1) * vxy(1) + uxy(2) * vxy(2) + uyy(2) * vyy(2)) + mu * iota^2 * ((uxy(1) + uxx(2)) * (vxy(1) + vxx(2)) + (uyy(1) + uxy(2)) * (vyy(1) + vxy(2)));
x = zeros(Nt*18^2, 1);
y = zeros(Nt*18^2, 1);
v = zeros(Nt*18^2, 1);
idx = @(k, i) (k - 1) * 18^2 + (i - 1) * 18 + (1: 18);
for k = 1: Nt
    K = zeros(18);
    v1 = vertices(mesh(k, 1), :);
    v2 = vertices(mesh(k, 3), :);
    v3 = vertices(mesh(k, 5), :);
    sgn1 = sign(mesh(k, 5) - mesh(k, 3));
    sgn2 = sign(mesh(k, 1) - mesh(k, 5));
    sgn3 = sign(mesh(k, 3) - mesh(k, 1));
    jacobi = T(v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
    for l = 1: size(wg, 1)
        vx = phix(ld(l, 1), ld(l, 2), ld(l, 3), sgn1, sgn2, sgn3, v1(1), v2(1), v3(1), v1(2), v2(2), v3(2)); %[18, 2]
        vy = phiy(ld(l, 1), ld(l, 2), ld(l, 3), sgn1, sgn2, sgn3, v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
        vxx = phixx(ld(l, 1), ld(l, 2), ld(l, 3), sgn1, sgn2, sgn3, v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
        vxy = phixy(ld(l, 1), ld(l, 2), ld(l, 3), sgn1, sgn2, sgn3, v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
        vyy = phiyy(ld(l, 1), ld(l, 2), ld(l, 3), sgn1, sgn2, sgn3, v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
        for j = 1: 18
            for i = 1: 18
                K(i, j) = K(i, j) + wg(l) * jacobi * a(vx(j, :), vy(j, :), vxx(j, :), vxy(j, :), vyy(j, :), vx(i, :), vy(i, :), vxx(i, :), vxy(i, :), vyy(i, :));
            end
        end
    end
    for i = 1: 18
        x(idx(k, i)) = mesh(k, i);
        y(idx(k, i)) = mesh(k, :);
        v(idx(k, i)) = K(i, :);
    end
end
A = sparse(x, y, v, Nv, Nv);
b = zeros(Nv, 1);
c = zeros(Nv, 1);
for k = 1: Nt
    v1 = vertices(mesh(k, 1), :);
    v2 = vertices(mesh(k, 3), :);
    v3 = vertices(mesh(k, 5), :);
    sgn1 = sign(mesh(k, 5) - mesh(k, 3));
    sgn2 = sign(mesh(k, 1) - mesh(k, 5));
    sgn3 = sign(mesh(k, 3) - mesh(k, 1));
    jacobi = T(v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
    for l = 1: size(wg, 1)
        x = v1(1) * ld(l, 1) + v2(1) * ld(l, 2) + v3(1) * ld(l, 3);
        y = v1(2) * ld(l, 1) + v2(2) * ld(l, 2) + v3(2) * ld(l, 3);
        v = phi(ld(l, 1), ld(l, 2), ld(l, 3), sgn1, sgn2, sgn3, v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
        b(mesh(k, 1: 18)) = b(mesh(k, 1: 18)) + wg(l) * jacobi * v * f(x, y);
    end
end
c(free) = A(free, free) \ b(free);
error = 0;
normu = 0;
for k = 1: Nt
    v1 = vertices(mesh(k, 1), :);
    v2 = vertices(mesh(k, 3), :);
    v3 = vertices(mesh(k, 5), :);
    sgn1 = sign(mesh(k, 5) - mesh(k, 3));
    sgn2 = sign(mesh(k, 1) - mesh(k, 5));
    sgn3 = sign(mesh(k, 3) - mesh(k, 1));    
    jacobi = T(v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
    for l = 1: size(wg, 1)
        vx = phix(ld(l, 1), ld(l, 2), ld(l, 3), sgn1, sgn2, sgn3, v1(1), v2(1), v3(1), v1(2), v2(2), v3(2)); %[18, 2]
        vy = phiy(ld(l, 1), ld(l, 2), ld(l, 3), sgn1, sgn2, sgn3, v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
        vxx = phixx(ld(l, 1), ld(l, 2), ld(l, 3), sgn1, sgn2, sgn3, v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
        vxy = phixy(ld(l, 1), ld(l, 2), ld(l, 3), sgn1, sgn2, sgn3, v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
        vyy = phiyy(ld(l, 1), ld(l, 2), ld(l, 3), sgn1, sgn2, sgn3, v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
        x = v1(1) * ld(l, 1) + v2(1) * ld(l, 2) + v3(1) * ld(l, 3);
        y = v1(2) * ld(l, 1) + v2(2) * ld(l, 2) + v3(2) * ld(l, 3);
        ex = ux(x, y) - vx' * c(mesh(k, 1: 18));
        ey = uy(x, y) - vy' * c(mesh(k, 1: 18));
        exx = uxx(x, y) - vxx' * c(mesh(k, 1: 18));
        exy = uxy(x, y) - vxy' * c(mesh(k, 1: 18));
        eyy = uyy(x, y) - vyy' * c(mesh(k, 1: 18));
        error = error + wg(l) * jacobi * a(ex, ey, exx, exy, eyy, ex, ey, exx, exy, eyy);
        normu = normu + wg(l) * jacobi * a(ux(x, y), uy(x, y), uxx(x, y), uxy(x, y), uyy(x, y), ux(x, y), uy(x, y), uxx(x, y), uxy(x, y), uyy(x, y));
    end
end
    error = (error / normu)^(1/2);
end