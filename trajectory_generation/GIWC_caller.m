clc 
clear all
close all
surfaces= get_surf('flat',3);
%surfaces= generate_box([0,0,0], 0.55, 0.4);
CoM=[0;0;10];
% figure 
% answ1=gravito_inertial_wrench_corretto(surfaces(1:2),CoM);
% axis equal
%figure

tic
for i = 1:100
    answ2=gravito_inertial_wrench(surfaces,CoM);
end
toc

%axis equal

goal_x = [0.4267,0.0964,0.2038,0.1633,-0.0045,-0.2123,-0.4232,-0.6161,-0.6965]