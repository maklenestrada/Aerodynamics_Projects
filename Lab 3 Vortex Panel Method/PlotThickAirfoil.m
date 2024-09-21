function PlotThickAirfoil(c,alpha,V_inf,p_inf,rho_inf,N,xx)
%Plot a thick airfoil
% Use the vortex panel method to plot a thick airfoil
%
% Author: Maklen A. Estrada
% Collaborators: Ian Wong, Mathew Pabin
% Date: 2/27/2022
%%
%Grid conditions
    Xmax = (3*c)/2;
    Xmin = -c/2;
    Ymax = c/2;
    Ymin = -c/2;
%Bound
    Xbound = [Xmin Xmax];
    Ybound = [Ymin Ymax];
    bound = [Xbound Ybound];
%Step size
    Xn = 100;
    Yn = 100;
%% NACA formula calculation 
%Thickness
    t = xx/100;
%Xcoordinates of airfoil
    XB_h1 = linspace(c,0,N/2);
    XB_h2 = flip(XB_h1(1:length(XB_h1)-1));
    XB = [XB_h1 XB_h2];
%Ycoordinates of airfoil 
    YB_h1 = (t/0.2)*c .* (  0.2969*sqrt(XB_h1/c) - 0.126*(XB_h1/c) - 0.3516*(XB_h1/c).^2 + 0.2843*(XB_h1/c).^3 - 0.1036*(XB_h1/c).^4); 
    YB_h2 = (t/0.2)*c .* (  0.2969*sqrt(XB_h2/c) - 0.126*(XB_h2/c) - 0.3516*(XB_h2/c).^2 + 0.2843*(XB_h2/c).^3 - 0.1036*(XB_h2/c).^4);
    YB = [-YB_h1 YB_h2];
%Grid
    [x,y] = meshgrid(linspace(Xmin,Xmax,Xn),linspace(Ymin,Ymax,Yn));
%Calling Vortex Panel Function
    [Cl,Cp,Circulation,X,Y] = Vortex_Panel(XB,YB,V_inf,rad2deg(alpha),0);
%% Phi and Psi Calculation
%VP for Vortex Flow and SL for Vortex flow
Phi_v = 0;
Psi_v = 0;
for i = 1:length(X)
    rad = sqrt((x-X(i)).^2 + (y-Y(i)).^2);
    theta = mod(atan2(y-Y(i),x-X(i)),2*pi);
    Phi_v = Phi_v - (Circulation(i)*theta)/(2*pi);
    Psi_v = Psi_v + (Circulation(i)*log(rad))/(2*pi);
end
% S.L (psi) for uniform flow
Psi_u = V_inf*(y*cos(alpha)-x*sin(alpha));
%Vortex and Uniform Flow
Phi = Phi_v + Phi_u;
% V.P (phi) for uniform flow
Phi_u = V_inf*(x*cos(alpha)-y*sin(alpha));
%Vortex and Uniform Flow
Psi = Psi_v + Psi_u;
%% Pressure
%Dynamic Pressure
q = (rho_inf*(V_inf^2))/2;
%Pressure Calculation 
Cp = 1-(gradient(Phi,(Xmax-Xmin)/Xn,(Ymax-Ymin)/Yn)./V_inf).^2;
Press = p_inf + Cp + q;
%% Plots 
%Streamline Plot
figure
contourf(x,y,Psi,90)
hold on 
%Airfoil
plot(XB,YB,'k--','linewidth',3)
fill(XB,YB,'w')
axis(bound)
ylabel('Y')
xlabel('X')
title(['Stream Lines for ' num2str(N) ' vorticies'])
hold off 

%Equipotential Plot
figure
contourf(x,y,Phi,75)
hold on 
%Airfoil
plot(XB,YB,'k--','linewidth',3)
fill(XB,YB,'w')
axis(bound)
ylabel('Y')
xlabel('X')
title(['Equipotential Lines for ' num2str(N) ' vorticies'])
hold off 

%Pressure Contour
figure
contourf(x,y,Press,80)
hold on 
%Airfoil
plot(XB,YB,'k--','linewidth',3)
hold on
fill(XB,YB,'w')
axis(bound)
ylabel('Y')
xlabel('X')
title(['Pressure Contour Lines for ' num2str(N) ' vorticies'])
hold off
end