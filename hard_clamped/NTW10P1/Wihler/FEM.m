function [erroru, errorp] = FEM(n, iota, lambda, mu, u, ux, uy, uxx, uxy, uyy, divu, phix, phiy, phixx, phixy, phiyy, p, px, py, T)
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
[vertices, edges, mesh] = init_mesh(n);
not_bdr = prod(vertices .* (1 - vertices), 2);
free = not_bdr ~= 0;
Nv = size(vertices, 1);
Nt = size(mesh, 1);
a = @(ux, uy, uxx, uxy, uyy, vx, vy, vxx, vxy, vyy) 2 * mu * (ux(1) * vx(1) + uy(2) * vy(2)) + mu * (uy(1) + ux(2)) * (vy(1) + vx(2)) + 2 * mu * iota^2 * (uxx(1) * vxx(1) + uxy(1) * vxy(1) + uxy(2) * vxy(2) + uyy(2) * vyy(2)) + mu * iota^2 * ((uxy(1) + uxx(2)) * (vxy(1) + vxx(2)) + (uyy(1) + uxy(2)) * (vyy(1) + vxy(2)));
b = @(vx, vy, vxx, vxy, vyy, q, qx, qy) (vx(1) + vy(2)) * q + iota^2 * ((vxx(1) + vxy(2)) * qx + (vxy(1) + vyy(2)) * qy);
C = @(p, px, py, q, qx, qy) p * q + iota^2 * (px * qx + py * qy);
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
                K(i+20, j+20) = K(i+20, j+20) - wg(l) * jacobi / lambda * C(q(j), qx(j), qy(j), q(i), qx(i), qy(i));
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
c = zeros(Nv, 1);
for k = 1: n
    v1 = vertices(edges(k, 1), :); %y=0
    v3 = vertices(edges(k, 2), :);
    sgn = sign(edges(k, 2) - edges(k, 1));
    c(edges(k, 1): edges(k, 1)+1) = u(v1(1), 0); 
    c(edges(k, 3): edges(k, 3)+1) = integral(@(x)u(x, 0), v1(1), v3(1), 'ArrayValued', true) / (v3(1) - v1(1));
    c(edges(k, 4): edges(k, 4)+1) = -sgn * integral(@(x)uy(x, 0), v1(1), v3(1), 'ArrayValued', true) / (v3(1) - v1(1));
    c(edges(k, 5)) = lambda *  divu(v1(1), 0);
    v1 = vertices(edges(n+k, 1), :); %x=1
    v3 = vertices(edges(n+k, 2), :);
    sgn = sign(edges(n+k, 2) - edges(n+k, 1));
    c(edges(n+k, 1): edges(n+k, 1)+1) = u((v1(2)^2+1)^(1/2), atan(v1(2))); 
    c(edges(n+k, 3): edges(n+k, 3)+1) = integral(@(y)u((y.^2+1).^(1/2), atan(y)), v1(2), v3(2), 'ArrayValued', true) / (v3(2) - v1(2));
    c(edges(n+k, 4): edges(n+k, 4)+1) = sgn * integral(@(y)ux((y.^2+1).^(1/2), atan(y)), v1(2), v3(2), 'ArrayValued', true) / (v3(2) - v1(2));
    c(edges(n+k, 5)) = lambda *  divu((v1(2)^2+1)^(1/2), atan(v1(2)));
    v1 = vertices(edges(2*n+k, 1), :); %y=1
    v3 = vertices(edges(2*n+k, 2), :);
    sgn = sign(edges(2*n+k, 2) - edges(2*n+k, 1));
    c(edges(2*n+k, 1): edges(2*n+k, 1)+1) = u((v1(1)^2+1)^(1/2), atan(1/v1(1))); 
    c(edges(2*n+k, 3): edges(2*n+k, 3)+1) = integral(@(x)u((x.^2+1).^(1/2), atan(1/x)), v1(1), v3(1), 'ArrayValued', true) / (v3(1) - v1(1));
    c(edges(2*n+k, 4): edges(2*n+k, 4)+1) = sgn * integral(@(x)uy((x.^2+1).^(1/2), atan(1/x)), v1(1), v3(1), 'ArrayValued', true) / (v3(1) - v1(1));
    c(edges(2*n+k, 5)) = lambda * divu((v1(1)^2+1)^(1/2), atan(1/v1(1)));
    v1 = vertices(edges(3*n+k, 1), :); %x=0
    v3 = vertices(edges(3*n+k, 2), :);
    sgn = sign(edges(3*n+k, 2) - edges(3*n+k, 1));
    c(edges(3*n+k, 1): edges(3*n+k, 1)+1) = u(v1(2), pi/2); 
    c(edges(3*n+k, 3): edges(3*n+k, 3)+1) = integral(@(y)u(y, pi/2), v1(2), v3(2), 'ArrayValued', true) / (v3(2) - v1(2));
    c(edges(3*n+k, 4): edges(3*n+k, 4)+1) = -sgn * integral(@(y)ux(y, pi/2), v1(2), v3(2), 'ArrayValued', true) / (v3(2) - v1(2));
    c(edges(3*n+k, 5)) = lambda * divu(v1(2), pi/2);
end
b = b - A * c;
c(free) = A(free, free) \ b(free);
erroru = 0;
normu = 0;
errorp = 0;
normp = 0;
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
    for l = 1: size(wg, 1)
        q = p(ld(l, 1), ld(l, 2), ld(l, 3));
        vx = phix(ld(l, 1), ld(l, 2), ld(l, 3), sgn1, sgn2, sgn3, v1(1), v2(1), v3(1), v1(2), v2(2), v3(2)); %[20, 2]
        vy = phiy(ld(l, 1), ld(l, 2), ld(l, 3), sgn1, sgn2, sgn3, v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
        vxx = phixx(ld(l, 1), ld(l, 2), ld(l, 3), sgn1, sgn2, sgn3, v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
        vxy = phixy(ld(l, 1), ld(l, 2), ld(l, 3), sgn1, sgn2, sgn3, v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
        vyy = phiyy(ld(l, 1), ld(l, 2), ld(l, 3), sgn1, sgn2, sgn3, v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
        x = v1(1) * ld(l, 1) + v2(1) * ld(l, 2) + v3(1) * ld(l, 3);
        y = v1(2) * ld(l, 1) + v2(2) * ld(l, 2) + v3(2) * ld(l, 3);
        rho = (x^2 + y^2)^(1/2);
        theta = atan(y/x);
        ex = ux(rho, theta) - vx' * c(mesh(k, 1: 20));
        ey = uy(rho, theta) - vy' * c(mesh(k, 1: 20));
        exx = uxx(rho, theta) - vxx' * c(mesh(k, 1: 20));
        exy = uxy(rho, theta) - vxy' * c(mesh(k, 1: 20));
        eyy = uyy(rho, theta) - vyy' * c(mesh(k, 1: 20));
        erroru = erroru + wg(l) * jacobi * a(ex, ey, exx, exy, eyy, ex, ey, exx, exy, eyy);
        normu = normu + wg(l) * jacobi * a(ux(rho, theta), uy(rho, theta), uxx(rho, theta), uxy(rho, theta), uyy(rho, theta), ux(rho, theta), uy(rho, theta), uxx(rho, theta), uxy(rho, theta), uyy(rho, theta));
        vx = lambda * ux(rho, theta);
        vy = lambda * uy(rho, theta);
        vxx = lambda * uxx(rho, theta);
        vxy = lambda * uxy(rho, theta);
        vyy = lambda * uyy(rho, theta);
        e = vx(1) + vy(2) - q' * c(mesh(k, 21: 23));
        ex = vxx(1) + vxy(2) - qx' * c(mesh(k, 21: 23));
        ey = vxy(1) + vyy(2) - qy' * c(mesh(k, 21: 23));
        errorp = errorp + wg(l) * jacobi * C(e, ex, ey, e, ex, ey);
        normp = normp + wg(l) * jacobi * C(vx(1)+vy(2), vxx(1)+vxy(2), vxy(1)+vyy(2), vx(1)+vy(2), vxx(1)+vxy(2), vxy(1)+vyy(2));
    end
end
    erroru = (erroru / normu)^(1/2);
    errorp = (errorp / normp)^(1/2);    
end