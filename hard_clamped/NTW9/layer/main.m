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
%[phi, phix, phiy, phixx, phixy, phiyy, T] = init_fespace;
load('basis_func');
error = zeros(1, 4);
rate = zeros(1, 4);
for i = 1: 2
    error(i) = FEM(2^(i+2), iota, lambda, mu, ux, uy, uxx, uxy, uyy, f, phi, phix, phiy, phixx, phixy, phiyy, T);
end
for i = 2: 2
    rate(i) = log2(error(i-1)/error(i));
end
fprintf("%.0e & %.3e & %.3e & %.3e & %.3e\\\\\n", iota, error(1), error(2), error(3), error(4));
fprintf("rate & & %.2f & %.2f & %.2f\\\\\n", rate(2), rate(3), rate(4));