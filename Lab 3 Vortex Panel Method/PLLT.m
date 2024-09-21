function [e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N)
%This function calculate the A coefficients using the Prandtls equation
%and calcualtes e, c_L, and c_Di
%
%Author: Maklen A. Estrada
%Collaborators: Ian Wong
%Date: November 6th, 2022
for i = 1:N
    theta(i) = (i*pi)/(2*N);
end
%Aspect Ratio
S = b*((c_t+c_r)/2);
AR = b^2/S;

%Chord calculation
chord = c_r + (c_t - c_r)*cos(theta);
%a_o calculation
a_o = a0_r + (a0_t - a0_r)*cos(theta);
%AoA L=0 calculation
alpha_0 = aero_r + (aero_t - aero_r)*cos(theta);
%Geometric AoA calculation
alpha_geo = geo_r + (geo_t - geo_r)*cos(theta);

%Solving for A coefficients
for i = 1:N
lhs(i) = alpha_geo(i) - alpha_0(i);%Left Hand Side of Prandtls eq
for j = 1:N
    n = (2*j - 1);
    rhs(i,j) = ((4*b*sin(n*theta(i))/(a_o(i)*chord(i))) + ((n*sin(n*theta(i)))/sin(theta(i))));%right hand side of P eq
end
end
A_coeff = inv(rhs)*lhs';

%Drag Factor 
A1 = A_coeff(1);
delta = 0;
for j = 2:N
    n = (2*j - 1);
    delta = delta + ((A_coeff(j)/A1)^2)*n;
end

%Lift Calculation 
c_L = A1*pi*AR;
%Span Efficiancy
e = (1+delta)^-1;
%Coeff of induced drag
c_Di = (c_L^2)/(pi*e*AR);
end