clear;
nu = 0.3;
iota = 1e-6;
E = 1;
lambda = E * nu / (1 + nu) / (1 - 2 * nu);
mu = E / 2 / (1 + nu);
syms x y;
u = [-sin(pi*x)^2*sin(2*pi*y); sin(2*pi*x)*sin(pi*y)^2];
divu = diff(u(1), x) + diff(u(2), y);
f = mu * (diff(diff(u, x), x) + diff(diff(u, y), y)) + (lambda + mu) * [diff(divu, x); diff(divu, y)];
f = -f;
ux = matlabFunction(simplify(diff(u, x)));
uy = matlabFunction(simplify(diff(u, y)));
uxx = matlabFunction(simplify(diff(diff(u, x), x)));
uxy = matlabFunction(simplify(diff(diff(u, x), y)));
uyy = matlabFunction(simplify(diff(diff(u, y), y)));
f = matlabFunction(simplify(f));
%[phi, phix, phiy, phixx, phixy, phiyy, p, px, py, T] = init_fespace;
load('basis_func');
erroru = zeros(1, 4);
rateu = zeros(1, 4);
errorp = zeros(1, 4);
ratep = zeros(1, 4);
for i = 1: 2
    [erroru(i), errorp(i)] = FEM(2^(i+2), iota, lambda, mu, ux, uy, uxx, uxy, uyy, f, phi, phix, phiy, phixx, phixy, phiyy, p, px, py, T);
end
for i = 2: 2
    rateu(i) = log2(erroru(i-1)/erroru(i));
    ratep(i) = log2(errorp(i-1)/errorp(i));
end
fprintf("%.0e & %.3e & %.3e & %.3e & %.3e\\\\\n", iota, erroru(1), erroru(2), erroru(3), erroru(4));
fprintf("rate & & %.2f & %.2f & %.2f\\\\\n", rateu(2), rateu(3), rateu(4));
%fprintf("%.0e & %.3e & %.3e & %.3e & %.3e\\\\\n", iota, errorp(1), errorp(2), errorp(3), errorp(4));
%fprintf("rate & & %.2f & %.2f & %.2f\\\\\n", ratep(2), ratep(3), ratep(4));