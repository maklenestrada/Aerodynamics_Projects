%% ASEN 3111 - Computational Assignment 1 - Main
%Maklen A. Estrada 
%Septeber 18 2021
clear all
clc
%% Problem 1 

%Coefficient of pressure function 
syms theta
Cp(theta) = (5/4) - 4*sin(theta)^2 - 2*sin(theta);

%Coef of lift and Coef of drag functions
Cl_fun(theta) = (-1/2)*Cp*sin(theta);
Cd_fun(theta) = (-1/2)*Cp*cos(theta);

%Analiticly integrting Cl and Cd
Cl = int(Cl_fun,[0 2*pi]);
Cd = int(Cd_fun,[0 2*pi]);

fprintf('Probelm 1')
fprintf('\n')
fprintf('Sectional Coefficient of Lift: %f',Cl)
fprintf('\n')
fprintf('Sectional Drag Coefficient: %f', Cd)
fprintf('\n')

%% Trapezoid Rule 
N = 5000;
x = linspace(0,2*pi,N);

%Coefficient of Lift calculation 
Cl_Trap = zeros(1,N);
for i = 1:N-1
Cl_Trap(i) = ( x(i+1) - x(i) ) * (...
    sin(x(i+1))^3 + sin(x(i))^3 +...
    (1/2)*( sin(x(i+1))^2 + sin(x(i))^2) + ...
    (-3/16)*( sin(x(i+1) + sin(x(i)))));
end

for i = 1:length(Cl_Trap)
Cl_T = sum(Cl_Trap(1:i));
Cl_TS(i) = sum(Cl_Trap(1:i));
end

%Coeficeint of drag calculation
Cd_Trap = zeros(1,N);
for i = 1:N-1
Cd_Trap(i) = ( x(i+1) - x(i)) * (...
    sin(x(i+1))^2*cos(x(i+1)) + sin(x(i))^2*cos(x(i)) +...
    (1/2)*( sin(x(i+1))*cos(x(i+1)) + sin(x(i))*cos(x(i))) +...
    (-3/16)*( cos(x(i+1)) + cos(x(i))) );
end

for i = 1:length(Cd_Trap)
Cd_T = sum(Cd_Trap(1:i));
Cd_TS(i) = sum(Cd_Trap(1:i));
end

figure 
plot(1:N,Cl_Trap,'k','LineWidth',2)
xlabel('Number of Panels')
ylabel('Coefficient of Lift')
title('Trapezoid Rule')
grid on

figure
plot(1:N,Cd_Trap,'b','LineWidth',2)
xlabel('Number of Panels')
ylabel('Coefficient of Drag')
title('Trapezoid Rule')
grid on

%% Simpsons Rule
N = 5000;
x = linspace(0,2*pi,N);

%Coefficient of Lift calculation 
for i = 2:(N-1)
Cl_Simp(i) = (pi/(3*N)) * (...
    4*sin(x(i-1))^3 + sin(x(i-1))^2 + (-3/8)*sin(x(i-1)) + ...
    4*(4*sin(x(i))^3 + sin(x(i))^2 + (-3/8)*sin(x(i))) +...
    4*sin(x(i+1))^3 + sin(x(i+1))^2 + (-3/8)*sin(x(i+1)));
end

for i = 1:length(Cl_Simp)
Cl_S = sum(Cl_Simp(1:i));
Cl_SS(i) = sum(Cl_Simp(1:i));
end

%Coeficeint of drag calculation
for i = 2:(N-1)
Cd_Simp(i) = (pi/(3*N)) * (...
    4*sin(x(i-1))^2*cos(x(i-1)) + sin(x(i-1))*cos(x(i-1)) + (-3/8)*cos(x(i-1)) +...
    4*(4*sin(x(i))^2*cos(x(i)) + sin(x(i))*cos(x(i)) + (-3/8)*cos(x(i))) +...
    4*sin(x(i+1))^2*cos(x(i+1)) + sin(x(i+1))*cos(x(i+1)) + (-3/8)*cos(x(i+1)));
end

for i = 1:length(Cd_Simp)
Cd_S = sum(Cd_Simp(1:i));
Cd_SS(i) = sum(Cd_Simp(1:i));
end

figure 
plot(1:(N-1),Cl_Simp,'k','LineWidth',2)
xlabel('Number of Panels')
ylabel('Coefficient of Lift')
title('Simpson Rule')
grid on

figure
plot(1:(N-1),Cd_Simp,'b','LineWidth',2)
xlabel('Number of Panels')
ylabel('Coefficient of Drag')
title('Simpson Rule')
grid on

%% Trapezoid Error 
Cl_Trap_Err =  abs(Cl_Trap - Cl)./Cl_Trap;
Cl_Tpan = find(Cl_Trap_Err < .01);
fprintf('Number of panels needed for sectional lift coefficient with 1 percent relative error using Trapezoidal Rule: %0.0f \n', Cl_Tpan(end))
%% Simpson Error 
Cl_Simp_Err =  abs(Cl_Simp - Cl)./Cl_Simp;
Cl_Span = find(Cl_Simp_Err < .01);
fprintf('Number of panels needed for sectional lift coefficient with 1 percent relative error using Simpson Rule: %0.0f \n', Cl_Span(end))
%% Problem 2
% Load Data
Data = load('Cp.mat');
Cp_U = Data.Cp_upper;
Cp_L = Data.Cp_lower;

%Constants
a = deg2rad(15);
c = 2;
V_inf = 30;
rho_inf = 1.225;
P_inf = 101.3e3;
t = 12/100;
N = 500000;
q_inf = (1/2)*rho_inf*(V_inf^2);

x = linspace(0,c,N);
y = ((t*c)/0.2) * (...
    0.2969*sqrt(x/c) - 0.1260*(x/c) - 0.3516*((x/c).^2) +...
    (0.2843*(x/c).^3) - 0.1036*((x/c).^4) ) ;

%Upper Surface Pressure
Upp = fnval(Cp_U,x/c);
P_upp = (Upp*q_inf) + P_inf;

%Lower Surface Pressure
Low = fnval(Cp_L,x/c);
P_low = (Low*q_inf) + P_inf; 

%Normal Force Calculation 
B_panel = (x(end) - x(1))/N;
Nor_Force = 0;
for i = 1:(N-1)
   N_Force(i) = Nor_Force;
   Nor_Force = Nor_Force +  B_panel * (...
       (P_low(i+1) + P_low(i))/2 -...
       (P_upp(i+1) + P_upp(i))/2 );
end

%Axial Force Calculation 
Ax_Force = 0;
for i = 1:(N-1)
    dT = y(i+1) - y(i);
    A_Force(i) = Ax_Force;
    Ax_Force = Ax_Force + dT* (...
        (P_low(i+1) + P_low(i))/2 +...
        (P_upp(i+1) + P_upp(i))/2 );
end

% % Lift and Drag Calculation
Lift = zeros(1,N);
Drag = zeros(1,N);
for i = 1:(N-1)
    Lift(i) = N_Force(i)*cos(a) - A_Force(i)*sin(a);
    Drag(i) = N_Force(i)*sin(a) + A_Force(i)*cos(a);
end

format longg
L = (Lift(end-1));
D = Drag(end-1);

fprintf('Probelm 2')
fprintf('\n')
fprintf('Lift: %f',L)
fprintf('\n')
fprintf('Drag: %f', D)
fprintf('\n')

%Error 
Lift_Err = abs(L - Lift)./ L;
L_E = find(Lift_Err < .01);
Drag_Err = abs(D - Drag)./ D;
D_E = find(Lift_Err < .01);

fprintf('Number of inegration points needed for Lift solution with less than 1 percent relative error: %0.0f \n', L_E(1))
fprintf('Number of inegration points needed for Drag solution with less than 1 percent relative error: %0.0f \n', D_E(1))
