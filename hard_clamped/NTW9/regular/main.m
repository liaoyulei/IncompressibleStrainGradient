clear;
nu = 0.3;
iota = 1;
E = 1;
lambda = E * nu / (1 + nu) / (1 - 2 * nu);
mu = E / 2 / (1 + nu);
syms x y;
u = [(cos(2*pi*x)-1)*(cos(4*pi*y)-1); (cos(4*pi*x)-1)*(cos(2*pi*y)-1)];
%u = [-sin(pi*x)^3*sin(2*pi*y)*sin(pi*y); sin(2*pi*x)*sin(pi*x)*sin(pi*y)^3];
divu = diff(u(1), x) + diff(u(2), y);
f = mu * (diff(diff(u, x), x) + diff(diff(u, y), y)) + (lambda + mu) * [diff(divu, x); diff(divu, y)];
f = iota^2 * (diff(diff(f, x), x) + diff(diff(f, y), y)) - f;
ux = matlabFunction(diff(u, x));
uy = matlabFunction(diff(u, y));
uxx = matlabFunction(diff(diff(u, x), x));
uxy = matlabFunction(diff(diff(u, x), y));
uyy = matlabFunction(diff(diff(u, y), y));
f = matlabFunction(f);
[phi, phix, phiy, phixx, phixy, phiyy, T] = init_fespace;
error = zeros(1, 4);
rate = zeros(1, 4);
for i = 1: 3
    error(i) = FEM(2^(i+2), iota, lambda, mu, ux, uy, uxx, uxy, uyy, f, phi, phix, phiy, phixx, phixy, phiyy, T);
end
for i = 2: 3
    rate(i) = log2(error(i-1)/error(i));
end
fprintf("%.0e & %.3e & %.3e & %.3e & %.3e\\\\\n", iota, error(1), error(2), error(3), error(4));
fprintf("rate & & %.2f & %.2f & %.2f\\\\\n", rate(2), rate(3), rate(4));