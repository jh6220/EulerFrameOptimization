function [K] = BeamF(L,h,alpha)
E = 70000;
G = E/(2*(1+0.3));
Iz = 6.35*h^3/12;
Ixy = h*6.35^3/12;
A = 6.35*h;
a = max(6.35,h);
b = min(6.35,h);
J = a*b^3*(1/3-0.21*b/a*(1-b^4/(12*a^4)));

k = zeros(12,12);
k(1,1) = E*A/L;
k(2,2) = 12*E*Iz/L^3;
k(3,3) = 12*E*Ixy/L^3;
k(4,4) = G*J/L;
k(5,5) = 4*E*Ixy/L;
k(6,6) = 4*E*Iz/L;
k(7,7) = E*A/L;
k(8,8) = 12*E*Iz/L^3;
k(9,9) = 12*E*Ixy/L^3;
k(10,10) = G*J/L;
k(11,11) = 4*E*Ixy/L;
k(12,12) = 4*E*Iz/L;
k(1,7) = -E*A/L; k(7,1) = -E*A/L;
k(2,6) = 6*E*Iz/L^2; k(6,2) = 6*E*Iz/L^2;
k(2,8) = -12*E*Iz/L^3; k(8,2) = -12*E*Iz/L^3;
k(2,12) = 6*E*Iz/L^2; k(12,2) = 6*E*Iz/L^2;
k(3,5) = -6*E*Ixy/L^2; k(5,3) = -6*E*Ixy/L^2;
k(3,9) = -12*E*Ixy/L^3; k(9,3) = -12*E*Ixy/L^3;
k(3,11) = -6*E*Ixy/L^2; k(11,3) = -6*E*Ixy/L^2;
k(4,10) = -G*J/L; k(10,4) = -G*J/L;
k(5,9) = 6*E*Ixy/L^2; k(9,5) = 6*E*Ixy/L^2;
k(5,11) = 2*E*Ixy/L; k(11,5) = 2*E*Ixy/L;
k(6,8) = -6*E*Iz/L^2; k(8,6) = -6*E*Iz/L^2;
k(6,12) = 2*E*Iz/L; k(12,6) = 2*E*Iz/L;
k(8,12) = -6*E*Iz/L^2; k(12,8) = -6*E*Iz/L^2;
k(9,11) = 6*E*Ixy/L^2; k(11,9) = 6*E*Ixy/L^2;

c = cos(alpha);
s = sin(alpha);

T = [c, s,  0,  0,  0,  0;
    -s, c,  0,  0,  0,  0;
    0,  0,  1,  0,  0,  0;
    0,  0,  0,  c,  s,  0;
    0,  0,  0,  -s, c,  0;
    0,  0,  0,  0,  0,  1];

T = [T,zeros(6,6);zeros(6,6),T];
K = transpose(T)*k*T;
end