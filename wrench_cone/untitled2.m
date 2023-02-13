clc
clear all

p11=[4, 2,   2]';
p21=[6, 2,   3]';
p31=[6, 6,   3]';
p41=[5, 8, 2.5]';
p51=[3, 6, 1.5]';
f11=[0, 1, 1]';
f21=[0, 0, 0.8]';
f31=[0.2, 0, 0.4]';
f41=[0, 0, 1]';
f51=[1, 0, 1]';
surface1.contact_pts=[p11 p21 p31 p41 p51]
surface1.contact_fs=[f11 f21 f31 f41 f51]
surface1.mu=0.65


p12=[-4, -8,  4]';
p22=[-6, -8,  4]';
p32=[-7, -4,  2]';
p42=[-5, -2,  1]';
p52=[-4, -4,  2]';
f12=[0, 0, 1]';
f22=[0, 0.1, 0.8]';
f32=[0, 0.2, 0.4]';
f42=[0, 0, 1]';
f52=[0.3, 0, 1]';
surface2.contact_pts=[p12 p22 p32 p42 p52]
surface2.contact_fs=[f12 f22 f32 f42 f52]
surface2.mu=0.65
surfaces=[surface1,surface2];

CoM=[0,0,10]'; %position of CoM
