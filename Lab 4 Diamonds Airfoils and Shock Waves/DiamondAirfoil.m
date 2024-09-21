function [c_l,c_dw] = DiamondAirfoil(M, alpha, epsilon1, epsilon2)
%This function takes in Mach, AoA, epsilon1, and epsilon. Then uses these
%paramaters to calculate the coefficient of lift an sectional wave-drag coefficient
%
% Author: Maklen Estrada
% Collaborators: Ian Wong
% Date: December 4th, 2022

%Constants 
gamma = 7/5;
type = 'Weak';

%Region 1
theta_1 = epsilon1-alpha;
beta_1 = ObliqueShockBeta(M,theta_1,gamma,type);
M_normal = M*sind(beta_1);
P_r1 = 1+((2*gamma)/(gamma+1))*(M_normal^2 - 1);
Mnormal_r1 = sqrt((1+((gamma-1)/2)*M_normal^2)/((gamma*M_normal^2)-((gamma-1)/2)));
M_r1 = Mnormal_r1/sind(beta_1-theta_1);


%Region 2
theta_2 = 2*epsilon2;
[~, V_r1, ~] = flowprandtlmeyer(gamma,abs(M_r1),'mach');
V_r2 = V_r1+theta_2;
[M_r2,~ ,~] = flowprandtlmeyer(gamma,V_r2,'nu');
Tratio_r2 = (1+(((gamma-1)/2)*M_r1^2))/(1+(((gamma-1)/2)*M_r2^2));
Pratio_r2 = Tratio_r2^(gamma/(gamma-1));
P_r2 = Pratio_r2 *P_r1;

%Region 3
theta_3 = epsilon1 + alpha;
beta_3 = ObliqueShockBeta(M,theta_3,gamma,type);
M_normal = M*sind(beta_3);
P_r3 = 1 + ((2*gamma)/(gamma+1))*(M_normal^2 - 1);
Mnormal_r3 = sqrt((1+((gamma - 1)/2)*M_normal^2)/((gamma-1)/2));
M_r3 = Mnormal_r3/sind(beta_3-theta_3);

%Region 4
theta4 = 2*epsilon2;
[~, V_r3, ~] = flowprandtlmeyer(gamma, M_r3, 'mach');
V_r4 = V_r3 +theta4;
[M4, ~, ~] = flowprandtlmeyer(gamma, V_r4, 'nu');
Tratio_r4 = (1+(((gamma-1)/2)*M_r3^2))/ (1+(((gamma-1)/2)*M4^2));
Pratio_r4 = Tratio_r4^(gamma/(gamma-1));
P_r4 = Pratio_r4 *P_r3;


%Calculating Coefficeint of lift and sectional wave-drag coefficeint
c_a = (1/(gamma*M^2))*(((P_r1+P_r3)*tand(epsilon1)*(2/3))+((-P_r2- P_r4)*tand(epsilon2)*(4/3)));
c_n = (1/(gamma*M^2))*((-P_r1+P_r3)*(2/3)+(-P_r2+P_r4)*(4/3));
c_l = c_n*cosd(alpha)-c_a*sind(alpha);
c_dw = c_n*sind(alpha)+c_a*cosd(alpha);
end