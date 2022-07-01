clear; %f=0
nu = 0.4999;
iota = 1e-6;
E = 1;
lambda = E * nu / (1 + nu) / (1 - 2 * nu);
mu = E / 2 / (1 + nu);
syms rho theta;
alpha = 1.5;
omega = 3/4 * pi;
c1 = -cos((alpha+1)*omega) / cos((alpha-1)*omega);
c2 = 2 * (lambda + 2 * mu) / (lambda + mu);
urho = rho^alpha / 2 / mu * (-(alpha + 1) * cos((alpha+1)*theta) + (c2 - alpha - 1) * c1 * cos((alpha-1)*theta));
utheta = rho^alpha / 2 / mu * ((alpha + 1) * sin((alpha+1)*theta) + (c2 + alpha - 1) * c1 * sin((alpha-1)*theta));
u = [cos(theta) * urho - sin(theta) * utheta; sin(theta) * urho + cos(theta) * utheta];
ux = diff(u, rho) * cos(theta) - diff(u, theta) * sin(theta) / rho;
uy = diff(u, rho) * sin(theta) + diff(u, theta) * cos(theta) / rho;
uxx = diff(ux, rho) * cos(theta) - diff(ux, theta) * sin(theta) / rho;
uxy = diff(ux, rho) * sin(theta) + diff(ux, theta) * cos(theta) / rho;
uyy = diff(uy, rho) * sin(theta) + diff(uy, theta) * cos(theta) / rho;
divu = simplify(ux(1) + uy(2));
u = matlabFunction(u); %rho theta
divu = matlabFunction(divu);
ux = matlabFunction(ux);
uy = matlabFunction(uy);
uxx = matlabFunction(uxx);
uxy = matlabFunction(uxy);
uyy = matlabFunction(uyy);
[phi, phix, phiy, phixx, phixy, phiyy, p, px, py, T] = init_fespace;
erroru = zeros(1, 4);
rateu = zeros(1, 4);
errorp = zeros(1, 4);
ratep = zeros(1, 4);
for i = 1: 4
    [erroru(i), errorp(i)] = FEM(2^(i+2), iota, lambda, mu, u, ux, uy, uxx, uxy, uyy, divu, phix, phiy, phixx, phixy, phiyy, p, px, py, T);
end
for i = 2: 4
    rateu(i) = log2(erroru(i-1)/erroru(i));
end
fprintf("%.0e & %.3e & %.3e & %.3e & %.3e\\\\\n", iota, erroru(1), erroru(2), erroru(3), erroru(4));
fprintf("rate & & %.2f & %.2f & %.2f\\\\\n", rateu(2), rateu(3), rateu(4));
%fprintf("%.0e & %.3e & %.3e & %.3e & %.3e\\\\\n", iota, errorp(1), errorp(2), errorp(3), errorp(4));
%fprintf("rate & & %.2f & %.2f & %.2f\\\\\n", ratep(2), ratep(3), ratep(4));