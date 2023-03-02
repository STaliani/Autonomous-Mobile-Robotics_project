clc 
clear all
close all
surfaces= get_surf1('four_points_nonflat');
CoM=[0;0;10];
% figure 
% answ1=gravito_inertial_wrench_corretto(surfaces(1:2),CoM);
% axis equal
figure
answ2=gravito_inertial_wrench(surfaces(1:2),CoM);
axis equal
