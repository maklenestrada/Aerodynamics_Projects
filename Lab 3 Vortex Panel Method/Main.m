%% ASEN 3111 - Computational Assignment 03 - Main
%This file is the main and calls the PLLT,Vortex_Panel NACA, and AirFoil
%functions
%
%Author: Maklen A. Estrada
%Collaborators: Ian Wong
%Date: November 6th, 2022

close all
clear all
clc
%% Problem 1 
c = 1;
xx = 12;
N = 100;
Vinf = 100;
alpha = deg2rad(10);
[Xb,Yb] = AirFoil(c,xx,N);
CL = Vortex_Panel(Xb,Yb,Vinf,alpha);

%Value of N = 1,000
CL_act = 0.0210;
CL_err = 1;
while CL_err > .01
    N = N+1;
    [Xb,Yb] = AirFoil(c,xx,N);
    CL = Vortex_Panel(Xb,Yb,Vinf,alpha);
    CL_err = abs((CL-CL_act)/CL_act);
end
fprintf('Problem 1')
fprintf('\n')
fprintf('The number of panels to get 1 %% error is %f ', N)
fprintf('\n')
%% Problem 2
V = 100;
alpha = linspace(-8,12,20);
% NACA0006
m0 = 0/100;
p0 = 0/10;
t0 = 6/100;
[x0,y0] = NACA(m0,p0,t0,c,N);%Getting Data from NACA Airfoils function
for i = 1:length(alpha)
[Cl0(i)] = Vortex_Panel(x0,y0,V,alpha(i));
end

% NACA 0012
m2 = 0/100;
p2 = 0/10;
t2 = 12/100;
[x2,y2] = NACA(m2,p2,t2,c,N);
for i = 1:length(alpha)
[Cl2(i)] = Vortex_Panel(x2,y2,V,alpha(i));
end

% NACA 0024 
m3 = 0/100;
p3 = 0/10;
t3 = 24/100;
[x3,y3] = NACA(m3,p3,t3,c,N);
for i = 1:length(alpha)
[Cl3(i)] = Vortex_Panel(x3,y3,V,alpha(i));
end

L1 = (Cl0(6) - Cl0(2))/(alpha(6) - alpha(2));
L2 = (Cl2(6) - Cl2(2))/(alpha(6) - alpha(2));
L3 = (Cl3(6) - Cl3(2))/(alpha(6) - alpha(2));
id1 = find(Cl0 == 0);
id2 = find(Cl2 == 0);
id3 = find(Cl3 == 0);

X = [L1 0;L2 0;L3 0];
X = array2table(X,'VariableNames',{'Lift Slope (1/rad)','Zero Lift AoA (deg)'},'RowNames',{'NACA 0006','NACA 0012','NACA 0024'});

fprintf('Problem 2')
fprintf('\n')
disp(X)
fprintf('\n')
figure(1)
plot(alpha,Cl0,'r')
hold on
plot(alpha,Cl2,'g')
hold on
plot(alpha,Cl3,'b')

title('Problem 2: \alpha vs. C_L')
legend('NACA 0006','NACA 0012','NACA 0024')
xlabel('\alpha (deg)')
ylabel('C_L')

%% Problem 3
V = 100;
alpha = linspace(-8,12,20);
% NACA 0012
m0 = 0/100;
p0 = 0/10;
t0 = 12/100;
[x0,y0] = NACA(m0,p0,t0,c,N);%Getting Data from NACA Airfoils function
for i = 1:length(alpha)
[Cl0(i)] = Vortex_Panel(x0,y0,V,alpha(i));
end

% NACA 2412
m2 = 2/100;
p2 = 4/10;
t2 = 12/100;
[x2,y2] = NACA(m2,p2,t2,c,N);
for i = 1:length(alpha)
[Cl2(i)] = Vortex_Panel(x2,y2,V,alpha(i));
end

% NACA 4412 
m3 = 4/100;
p3 = 4/10;
t3 = 12/100;
[x3,y3] = NACA(m3,p3,t3,c,N);
for i = 1:length(alpha)
[Cl3(i)] = Vortex_Panel(x3,y3,V,alpha(i));
end

L1 = (Cl0(6) - Cl0(2))/(alpha(6) - alpha(2));
L2 = (Cl2(6) - Cl2(2))/(alpha(6) - alpha(2));
L3 = (Cl3(6) - Cl3(2))/(alpha(6) - alpha(2));
id1 = find(Cl0 == 0);
id2 = find(Cl2 == 0);
id3 = find(Cl3 == 0);

X = [L1 alpha(9);L2 alpha(7);L3 alpha(5)];
X = array2table(X,'VariableNames',{'Lift Slope (1/rad)','Zero Lift AoA (deg)'},'RowNames',{'NACA 0012','NACA 2412','NACA 4424'});
fprintf('Problem 3')
fprintf('\n')
disp(X)
fprintf('\n')
figure(2)
plot(alpha,Cl0,'r')
hold on
plot(alpha,Cl2,'g')
hold on
plot(alpha,Cl3,'b')

title('Problem 3: \alpha vs. C_L')
legend('NACA 0006','NACA 0012','NACA 0024')
xlabel('\alpha (deg)')
ylabel('C_L')

%% Problem 5
%Constants 
V = 134.4;%ft/s
rho = .001756;
b = 33 + (4/12);
a0_t = (Cl0(6)-Cl0(2))/(alpha(6) - alpha(2));
a0_r = (Cl2(6)-Cl2(2))/(alpha(6) - alpha(2));
c_t = 3.7083;
c_r = 5.3333;
geo_t = 0;
geo_r = deg2rad(1);

for i=1:length(Cl2)
        aero_r = alpha(i);
    if Cl2(i) <=0
        break;
    end
end

for i=1:length(Cl0)
    aero_t = alpha(i);
    if Cl0(i) <=0
        break;
    end
end

[e, c_L, c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N);

S = b*((c_t+c_r)/2);
q = (1/2)*rho*V^2;
Lift = c_L*q*S;
Drag = c_Di*q*S;

fprintf('Problem 5')
fprintf('\n')
fprintf('Lift: %f (lbf)',Lift)
fprintf('\n')
fprintf('Induced Drag: %f (lbf)',Drag)
fprintf('\n')

%Lift error calculations 
cL = zeros(1,N);
for i = 1:N
    [~,cL(i),~] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,i);
end
cLerror = (c_L-cL);

%Error 10%
e1 = c_L*0.1;
numb1 = find(cLerror < e1);
n1 = numb1(1);
%Error 1%
e2 = c_L*0.01;
numb2 = find(cLerror < e2);
n2 = numb2(1);
%Error 1/10%
e3 = c_L*0.001;
numb3 = find(cLerror < e3);
n3 = numb3(1);

fprintf('Number of odd terms to obtain Lift solution within 10 %% relatibe error %f ',n1)
fprintf('\n')
fprintf('Number of odd terms to obtain Lift solution within 1 %% relatibe error %f ',n2)
fprintf('\n')
fprintf('Number of odd terms to obtain Lift solution within 0.1 %% relatibe error %f ',n3)
fprintf('\n')

%Induced Drag error calculations
cDi = zeros(1,N);
for i = 1:N
    [~,~,cDi(i)] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,i);
end
cDierror = (c_Di-cDi);

%Error 10%
e1 = c_Di*0.1;
numb1 = find(cDierror < e1);
n1 = numb1(1);
%Error 1%
e2 = c_Di*0.01;
numb2 = find(cDierror < e2);
n2 = numb2(1);
%Error 1/10%
e3 = c_Di*0.001;
numb3 = find(cDierror < e3);
n3 = numb3(1);

fprintf('Number of odd terms to obtain Induced Drag solution within 10 %% relatibe error %f ',n1)
fprintf('\n')
fprintf('Number of odd terms to obtain Induced Drag solution within 1 %% relatibe error %f ',n2)
fprintf('\n')
fprintf('Number of odd terms to obtain Induced Drag solution within 0.1 %% relatibe error %f ',n3)
