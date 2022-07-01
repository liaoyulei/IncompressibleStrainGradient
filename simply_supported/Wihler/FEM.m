function erroru = FEM(n, iota, lambda, mu, u, ux, uy, uxx, uxy, uyy, phix, phiy, phixx, phixy, phiyy, p, px, py, T)
wg = [0.144315607677787; 0.103217370534718; 0.103217370534718; 0.103217370534718; 0.032458497623198; 0.032458497623198; 0.032458497623198; 0.095091634267285; 0.095091634267285; 0.095091634267285; 0.027230314174435; 0.027230314174435; 0.027230314174435; 0.027230314174435; 0.027230314174435; 0.027230314174435];
ld = [
    1/3, 1/3, 1/3; 
    0.170569307751760, 0.170569307751760, 0.658861384496480;
    0.170569307751760, 0.658861384496480, 0.170569307751760;
    0.658861384496480, 0.170569307751760, 0.170569307751760;
    0.050547228317031, 0.050547228317031, 0.898905543365938;
    0.050547228317031, 0.898905543365938, 0.050547228317031;
    0.898905543365938, 0.050547228317031, 0.050547228317031;
    0.459292588292723, 0.459292588292723, 0.081414823414554;
    0.459292588292723, 0.081414823414554, 0.459292588292723;
    0.081414823414554, 0.459292588292723, 0.459292588292723;
    0.263112829634638, 0.728492392955404, 0.008394777409958;
    0.263112829634638, 0.008394777409958, 0.728492392955404;
    0.728492392955404, 0.263112829634638, 0.008394777409958;
    0.728492392955404, 0.008394777409958, 0.263112829634638;
    0.008394777409958, 0.728492392955404, 0.263112829634638;
    0.008394777409958, 0.263112829634638, 0.728492392955404
];
wg1 = [1/2; 1/2];
ld1 = ([-1/sqrt(3); 1/sqrt(3)] + 1) / 2;
[vertices, mesh] = init_mesh(n);
not_bdr = prod(vertices .* (1 - vertices), 2);
free = not_bdr ~= 0;
Nv = size(vertices, 1);
Nt = size(mesh, 1);
a = @(ux, uy, uxx, uxy, uyy, vx, vy, vxx, vxy, vyy) 2 * mu * (ux(1) * vx(1) + uy(2) * vy(2)) + mu * (uy(1) + ux(2)) * (vy(1) + vx(2)) + 2 * mu * iota^2 * (uxx(1) * vxx(1) + uxy(1) * vxy(1) + uxy(2) * vxy(2) + uyy(2) * vyy(2)) + mu * iota^2 * ((uxy(1) + uxx(2)) * (vxy(1) + vxx(2)) + (uyy(1) + uxy(2)) * (vyy(1) + vxy(2)));
b = @(vx, vy, vxx, vxy, vyy, q, qx, qy) (vx(1) + vy(2)) * q + iota^2 * ((vxx(1) + vxy(2)) * qx + (vxy(1) + vyy(2)) * qy);
c = @(p, px, py, q, qx, qy) p * q + iota^2 * (px * qx + py * qy);
B = zeros(Nv, 1);
C = zeros(Nv, 1);
x = zeros(Nt*23^2, 1);
y = zeros(Nt*23^2, 1);
v = zeros(Nt*23^2, 1);
idx = @(k, i) (k - 1) * 23^2 + (i - 1) * 23 + (1: 23);
for k = 1: Nt
    K = zeros(23);
    v1 = vertices(mesh(k, 1), :);
    v2 = vertices(mesh(k, 3), :);
    v3 = vertices(mesh(k, 5), :);
    v4 = vertices(mesh(k, 7), :);   
    sgn1 = sign(mesh(k, 5) - mesh(k, 3));
    sgn2 = sign(mesh(k, 1) - mesh(k, 5));
    sgn3 = sign(mesh(k, 3) - mesh(k, 1));
    qx = px(v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
    qy = py(v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
    jacobi = T(v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
    for l = 1: size(wg, 1)
        vx = phix(ld(l, 1), ld(l, 2), ld(l, 3), sgn1, sgn2, sgn3, v1(1), v2(1), v3(1), v1(2), v2(2), v3(2)); %[20, 2]
        vy = phiy(ld(l, 1), ld(l, 2), ld(l, 3), sgn1, sgn2, sgn3, v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
        vxx = phixx(ld(l, 1), ld(l, 2), ld(l, 3), sgn1, sgn2, sgn3, v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
        vxy = phixy(ld(l, 1), ld(l, 2), ld(l, 3), sgn1, sgn2, sgn3, v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
        vyy = phiyy(ld(l, 1), ld(l, 2), ld(l, 3), sgn1, sgn2, sgn3, v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
        q = p(ld(l, 1), ld(l, 2), ld(l, 3));
        for j = 1: 20
            for i = 1: 20
                K(i, j) = K(i, j) + wg(l) * jacobi * a(vx(j, :), vy(j, :), vxx(j, :), vxy(j, :), vyy(j, :), vx(i, :), vy(i, :), vxx(i, :), vxy(i, :), vyy(i, :));
            end
            for i = 1: 3
                K(i+20, j) = K(i+20, j) + wg(l) * jacobi * b(vx(j, :), vy(j, :), vxx(j, :), vxy(j, :), vyy(j, :), q(i), qx(i), qy(i));
            end
        end
        for j = 1: 3
            for i = 1: 20
                K(i, j+20) = K(i, j+20) + wg(l) * jacobi * b(vx(i, :), vy(i, :), vxx(i, :), vxy(i, :), vyy(i, :), q(j), qx(j), qy(j));
            end
            for i = 1: 3
                K(i+20, j+20) = K(i+20, j+20) - wg(l) * jacobi / lambda * c(q(j), qx(j), qy(j), q(i), qx(i), qy(i));
            end
        end
    end
    for i = 1: 23
        x(idx(k, i)) = mesh(k, i);
        y(idx(k, i)) = mesh(k, :);
        v(idx(k, i)) = K(i, :);
    end
    if free(mesh(k, 7)) == 0
        for l = 1: 2
            vx = phix(0, ld1(l), ld1(3-l), sgn1, sgn2, sgn3, v1(1), v2(1), v3(1), v1(2), v2(2), v3(2)); %[20, 2]
            vy = phiy(0, ld1(l), ld1(3-l), sgn1, sgn2, sgn3, v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
            rho = ((v2(1) * ld1(l) + v3(1) * ld1(3-l))^2 + (v2(2) * ld1(l) + v3(2) * ld1(3-l))^2)^0.5;
            theta = atan((v2(2) * ld1(l) + v3(2) * ld1(3-l)) / (v2(1) * ld1(l) + v3(1) * ld1(3-l)));
            vxx = uxx(rho, theta);
            vxy = uxy(rho, theta);
            vyy = uyy(rho, theta);
            sigmann = mu * [2*vxx(1), vxy(2)+vyy(1); vxx(2)+vxy(1), 2*vyy(2)] + lambda * [vxx(1)+vxy(2), 0; 0, vxy(1)+vyy(2)];
            if v4(2) == 0 %y=0
                B(mesh(k, 1: 20)) = B(mesh(k, 1: 20)) - wg1(l) * iota^2 / n * vy * sigmann(:, 2);
                C(mesh(k, [3, 7, 4, 8])) = u([v2(1); v4(1)], 0);
            elseif v4(1) == 1 %x=1
                B(mesh(k, 1: 20)) = B(mesh(k, 1: 20)) + wg1(l) * iota^2 / n * vx * sigmann(:, 1);
                C(mesh(k, [3, 7, 4, 8])) = u(([v2(2); v4(2)].^2+1).^0.5, atan([v2(2); v4(2)]));
            elseif v4(2) == 1 %y=1
                B(mesh(k, 1: 20)) = B(mesh(k, 1: 20)) + wg1(l) * iota^2 / n * vy * sigmann(:, 2);
                C(mesh(k, [3, 7, 4, 8])) = u(([v2(1); v4(1)].^2+1).^0.5, atan(1./[v2(1); v4(1)]));
            elseif v4(1) == 0 %x=0
                B(mesh(k, 1: 20)) = B(mesh(k, 1: 20)) - wg1(l) * iota^2 / n * vx * sigmann(:, 1);
                C(mesh(k, [3, 7, 4, 8])) = u([v2(2); v4(2)], pi/2);
            end
        end
    end
end
A = sparse(x, y, v, Nv, Nv);
B = B - A * C;
C(free) = A(free, free) \ B(free);
erroru = 0;
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
        vx = phix(ld(l, 1), ld(l, 2), ld(l, 3), sgn1, sgn2, sgn3, v1(1), v2(1), v3(1), v1(2), v2(2), v3(2)); %[20, 2]
        vy = phiy(ld(l, 1), ld(l, 2), ld(l, 3), sgn1, sgn2, sgn3, v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
        vxx = phixx(ld(l, 1), ld(l, 2), ld(l, 3), sgn1, sgn2, sgn3, v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
        vxy = phixy(ld(l, 1), ld(l, 2), ld(l, 3), sgn1, sgn2, sgn3, v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
        vyy = phiyy(ld(l, 1), ld(l, 2), ld(l, 3), sgn1, sgn2, sgn3, v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
        x = v1(1) * ld(l, 1) + v2(1) * ld(l, 2) + v3(1) * ld(l, 3);
        y = v1(2) * ld(l, 1) + v2(2) * ld(l, 2) + v3(2) * ld(l, 3);
        rho = (x^2 + y^2)^(1/2);
        theta = atan(y/x);
        ex = ux(rho, theta) - vx' * C(mesh(k, 1: 20));
        ey = uy(rho, theta) - vy' * C(mesh(k, 1: 20));
        exx = uxx(rho, theta) - vxx' * C(mesh(k, 1: 20));
        exy = uxy(rho, theta) - vxy' * C(mesh(k, 1: 20));
        eyy = uyy(rho, theta) - vyy' * C(mesh(k, 1: 20));
        erroru = erroru + wg(l) * jacobi * a(ex, ey, exx, exy, eyy, ex, ey, exx, exy, eyy);
        normu = normu + wg(l) * jacobi * a(ux(rho, theta), uy(rho, theta), uxx(rho, theta), uxy(rho, theta), uyy(rho, theta), ux(rho, theta), uy(rho, theta), uxx(rho, theta), uxy(rho, theta), uyy(rho, theta));
    end
end
    erroru = (erroru / normu)^(1/2);
end