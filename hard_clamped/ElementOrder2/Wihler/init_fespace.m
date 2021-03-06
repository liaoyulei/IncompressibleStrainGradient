function  [phi, phix, phiy, phixx, phixy, phiyy, p, px, py, T] = init_fespace
syms lambda1 lambda2 lambda3 x1 x2 x3 y1 y2 y3 sgn1 sgn2 sgn3
xi1 = x2 - x3;
xi2 = x3 - x1;
xi3 = x1 - x2;
eta1 = y2 - y3;
eta2 = y3 - y1;
eta3 = y1 - y2;
inner11 = eta1^2 + xi1^2;
inner12 = eta1 * eta2 + xi1 * xi2;
inner13 = eta1 * eta3 + xi1 * xi3;
inner22 = eta2^2 + xi2^2;
inner23 = eta2 * eta3 + xi2 * xi3;
inner33 = eta3^2 + xi3^2;
T2 = det([1, x1, y1; 1, x2, y2; 1, x3, y3]);
bK = lambda1 * lambda2 * lambda3;
b1 = lambda2 * lambda3;
b2 = lambda1 * lambda3;
b3 = lambda1 * lambda2;
phi = [
    lambda1*(3*lambda1-2)+60*bK*b1+30*inner12/inner22*bK*b2*(4*lambda2-1)+30*inner13/inner33*bK*b3*(4*lambda3-1)+180*bK^2;
    lambda2*(3*lambda2-2)+60*bK*b2+30*inner12/inner11*bK*b1*(4*lambda1-1)+30*inner23/inner33*bK*b3*(4*lambda3-1)+180*bK^2;
    lambda3*(3*lambda3-2)+60*bK*b3+30*inner13/inner11*bK*b1*(4*lambda1-1)+30*inner23/inner22*bK*b2*(4*lambda2-1)+180*bK^2;
    6*b1+90*bK*b1-90*bK*b2-90*bK*b3-900*bK^2;
    6*b2+90*bK*b2-90*bK*b1-90*bK*b3-900*bK^2;
    6*b3+90*bK*b3-90*bK*b1-90*bK*b2-900*bK^2;
    sgn1*30*T2/sqrt(inner11)*bK*b1*(4*lambda1-1);
    sgn2*30*T2/sqrt(inner22)*bK*b2*(4*lambda2-1);
    sgn3*30*T2/sqrt(inner33)*bK*b3*(4*lambda3-1);
    2520*bK^2;
];
lambda1x = eta1 / T2;
lambda2x = eta2 / T2;
lambda3x = eta3 / T2;
lambda1y = -xi1 / T2;
lambda2y = -xi2 / T2;
lambda3y = -xi3 / T2;
phi2(1: 2: 20, 1) = phi;
phi2(2: 2: 20, 2) = phi;
phi = phi2;
phix = lambda1x * diff(phi, lambda1) + lambda2x * diff(phi, lambda2) + lambda3x * diff(phi, lambda3);
phiy = lambda1y * diff(phi, lambda1) + lambda2y * diff(phi, lambda2) + lambda3y * diff(phi, lambda3);
phixx = lambda1x * diff(phix, lambda1) + lambda2x * diff(phix, lambda2) + lambda3x * diff(phix, lambda3);
phixy = lambda1y * diff(phix, lambda1) + lambda2y * diff(phix, lambda2) + lambda3y * diff(phix, lambda3);
phiyy = lambda1y * diff(phiy, lambda1) + lambda2y * diff(phiy, lambda2) + lambda3y * diff(phiy, lambda3);
phi = matlabFunction(phi);
phix = matlabFunction(phix);
phiy = matlabFunction(phiy);
phixx = matlabFunction(phixx);
phixy = matlabFunction(phixy);
phiyy = matlabFunction(phiyy);
p = matlabFunction([lambda1; lambda2; lambda3]);
px = matlabFunction([lambda1x; lambda2x; lambda3x]);
py = matlabFunction([lambda1y; lambda2y; lambda3y]);
T = matlabFunction(T2/2);
end