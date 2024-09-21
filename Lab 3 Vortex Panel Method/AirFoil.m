function [XB,YB] = AirFoil(c,xx,N)
%This function calculates the x and y coordinates
%
%Author: Maklen A. Estrada
%Date: November 6th, 2022
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
end