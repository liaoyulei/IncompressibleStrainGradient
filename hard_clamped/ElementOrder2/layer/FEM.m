function [erroru, errorp] = FEM(n, iota, lambda, mu, ux, uy, uxx, uxy, uyy, f, phi, phix, phiy, phixx, phixy, phiyy, p, px, py, T)
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
[vertices, mesh] = init_mesh(n);
not_bdr = prod(vertices .* (1 - vertices), 2);
free = not_bdr ~= 0;
Nv = size(vertices, 1);
Nt = size(mesh, 1);
a = @(ux, uy, uxx, uxy, uyy, vx, vy, vxx, vxy, vyy) 2 * mu * (ux(1) * vx(1) + uy(2) * vy(2)) + mu * (uy(1) + ux(2)) * (vy(1) + vx(2)) + 2 * mu * iota^2 * (uxx(1) * vxx(1) + uxy(1) * vxy(1) + uxy(2) * vxy(2) + uyy(2) * vyy(2)) + mu * iota^2 * ((uxy(1) + uxx(2)) * (vxy(1) + vxx(2)) + (uyy(1) + uxy(2)) * (vyy(1) + vxy(2)));
b = @(vx, vy, vxx, vxy, vyy, q, qx, qy) (vx(1) + vy(2)) * q + iota^2 * ((vxx(1) + vxy(2)) * qx + (vxy(1) + vyy(2)) * qy);
c = @(p, px, py, q, qx, qy) p * q + iota^2 * (px * qx + py * qy);
x = zeros(Nt*23^2, 1);
y = zeros(Nt*23^2, 1);
v = zeros(Nt*23^2, 1);
idx = @(k, i) (k - 1) * 23^2 + (i - 1) * 23 + (1: 23);
for k = 1: Nt
    K = zeros(23);
    v1 = vertices(mesh(k, 1), :);
    v2 = vertices(mesh(k, 3), :);
    v3 = vertices(mesh(k, 5), :);   
    sgn1 = sign(mesh(k, 5) - mesh(k, 3));
    sgn2 = sign(mesh(k, 1) - mesh(k, 5));
    sgn3 = sign(mesh(k, 3) - mesh(k, 1));
    qx = px(v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
    qy = py(v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
    jacobi = T(v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
    for l = 1: 16
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
end
A = sparse(x, y, v, Nv, Nv);
b = zeros(Nv, 1);
u = zeros(Nv, 1);
for k = 1: Nt
    v1 = vertices(mesh(k, 1), :);
    v2 = vertices(mesh(k, 3), :);
    v3 = vertices(mesh(k, 5), :);
    sgn1 = sign(mesh(k, 5) - mesh(k, 3));
    sgn2 = sign(mesh(k, 1) - mesh(k, 5));
    sgn3 = sign(mesh(k, 3) - mesh(k, 1));    
    jacobi = T(v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
    for l = 1: 16 
        v = phi(ld(l, 1), ld(l, 2), ld(l, 3), sgn1, sgn2, sgn3, v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
        x = v1(1) * ld(l, 1) + v2(1) * ld(l, 2) + v3(1) * ld(l, 3);
        y = v1(2) * ld(l, 1) + v2(2) * ld(l, 2) + v3(2) * ld(l, 3);
        b(mesh(k, 1: 20)) = b(mesh(k, 1: 20)) + wg(l) * jacobi * v * f(x, y);
    end
end
u(free) = A(free, free) \ b(free);
erroru = 0;
normu = 0;
errorp = 0;
for k = 1: Nt
    v1 = vertices(mesh(k, 1), :);
    v2 = vertices(mesh(k, 3), :);
    v3 = vertices(mesh(k, 5), :);
    sgn1 = sign(mesh(k, 5) - mesh(k, 3));
    sgn2 = sign(mesh(k, 1) - mesh(k, 5));
    sgn3 = sign(mesh(k, 3) - mesh(k, 1));
    jacobi = T(v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
    qx = px(v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
    qy = py(v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
    for l = 1: 16
        q = p(ld(l, 1), ld(l, 2), ld(l, 3));
        vx = phix(ld(l, 1), ld(l, 2), ld(l, 3), sgn1, sgn2, sgn3, v1(1), v2(1), v3(1), v1(2), v2(2), v3(2)); %[20, 2]
        vy = phiy(ld(l, 1), ld(l, 2), ld(l, 3), sgn1, sgn2, sgn3, v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
        vxx = phixx(ld(l, 1), ld(l, 2), ld(l, 3), sgn1, sgn2, sgn3, v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
        vxy = phixy(ld(l, 1), ld(l, 2), ld(l, 3), sgn1, sgn2, sgn3, v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
        vyy = phiyy(ld(l, 1), ld(l, 2), ld(l, 3), sgn1, sgn2, sgn3, v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
        x = v1(1) * ld(l, 1) + v2(1) * ld(l, 2) + v3(1) * ld(l, 3);
        y = v1(2) * ld(l, 1) + v2(2) * ld(l, 2) + v3(2) * ld(l, 3);
        ex = ux(x, y) - vx' * u(mesh(k, 1: 20));
        ey = uy(x, y) - vy' * u(mesh(k, 1: 20));
        exx = uxx(x, y) - vxx' * u(mesh(k, 1: 20));
        exy = uxy(x, y) - vxy' * u(mesh(k, 1: 20));
        eyy = uyy(x, y) - vyy' * u(mesh(k, 1: 20));
        erroru = erroru + wg(l) * jacobi * a(ex, ey, exx, exy, eyy, ex, ey, exx, exy, eyy);
        normu = normu + wg(l) * jacobi * a(ux(x, y), uy(x, y), uxx(x, y), uxy(x, y), uyy(x, y), ux(x, y), uy(x, y), uxx(x, y), uxy(x, y), uyy(x, y));
        vx = lambda * ux(x, y);
        vy = lambda * uy(x, y);
        vxx = lambda * uxx(x, y);
        vxy = lambda * uxy(x, y);
        vyy = lambda * uyy(x, y);
        e = vx(1) + vy(2) - q' * u(mesh(k, 21: 23));
        ex = vxx(1) + vxy(2) - qx' * u(mesh(k, 21: 23));
        ey = vxy(1) + vyy(2) - qy' * u(mesh(k, 21: 23));
        errorp = errorp + wg(l) * jacobi * c(e, ex, ey, e, ex, ey);
    end
end
erroru = (erroru / normu)^(1/2);
errorp = errorp^(1/2);
plotx = 0: 0.01: 1;
ploty = zeros(1, 101);
for k = 1: Nt
    v1 = vertices(mesh(k, 1), :);
    v2 = vertices(mesh(k, 3), :);
    v3 = vertices(mesh(k, 5), :);
    sgn1 = sign(mesh(k, 5) - mesh(k, 3));
    sgn2 = sign(mesh(k, 1) - mesh(k, 5));
    sgn3 = sign(mesh(k, 3) - mesh(k, 1));
    if (v1(2)<=0.5||v2(2)<=0.5||v3(2)<=0.5) && (v1(2)>=0.5||v2(2)>=0.5|| v3(2)>=0.5)
        for i = ceil(min([v1(1), v2(1), v3(1)]) / 0.01) : floor(max([v1(1), v2(1), v3(1)]) / 0.01)
            x = 0.01 * i;
            lambda1 = det([1, x, 0.5; 1, v2(1), v2(2); 1, v3(1), v3(2)]) / 2 / jacobi;
            lambda2 = det([1, v1(1), v1(2); 1, x, 0.5; 1, v3(1), v3(2)]) / 2 / jacobi;
            lambda3 = det([1, v1(1), v1(2); 1, v2(1), v2(2); 1, x, 0.5]) / 2 / jacobi;
            if lambda1 >=0 && lambda2 >= 0 && lambda3 >= 0
                vx = phi(lambda1, lambda2, lambda3, sgn1, sgn2, sgn3, v1(1), v2(1), v3(1), v1(2), v2(2), v3(2)); %[20, 2]
                y = vx' * u(mesh(k, 1: 20));
                ploty(i+1) = y(2);
            end
        end
    end
end
plot(plotx, ploty);
end