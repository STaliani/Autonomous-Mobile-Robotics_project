clc
clear all
addpath('casadi-windows-matlabR2016a-v3.5.5')
%addpath('C:/Users/anto1/Documents/MATLAB/casadi-windows-matlabR2016a-v3.5.5')
%addpath('C:/Users/saver/Documents/MATLAB/casadi-windows-matlabR2016a-v3.5.5')
import casadi.*

n_stance = 4;

p11=[-4, -8,   2]';
p21=[-6, -8,   3]';
p31=[-7, -4,   3.5]';
p41=[-5, -2,   2.5]';
p51=[-4, -4,   2]';
f11=[0, 1, 1]';
f21=[0, 0, 0.8]';
f31=[0.2, 0, 0.4]';
f41=[0, 0, 1]';
f51=[1, 0, 1]';
surface1.contact_pts=[p11 p21 p31 p41 p51];
surface1.contact_fs=[f11 f21 f31 f41 f51];
surface1.mu=0.65;
p12=[4, 2,  2]';
p22=[6, 2,  3]';
p32=[7, 6,  3.5]';
p42=[5, 8,  2.5]';
p52=[4, 6,  2]';
f12=[0, 0, 1]';
f22=[0, 0.1, 0.8]';
f32=[0, 0.2, 0.4]';
f42=[0, 0, 1]';
f52=[0.3, 0, 1]';
surface2.contact_pts=[p12 p22 p32 p42 p52];
surface2.contact_fs=[f12 f22 f32 f42 f52];
surface2.mu=0.65;


p13=[-4, 12,   2]';
p23=[-6, 12,   3]';
p33=[-7, 16,   3.5]';
p43=[-5, 18,   2.5]';
p53=[-4, 16,   2]';
f13=[0, 1, 1]';
f23=[0, 0, 0.8]';
f33=[0.2, 0, 0.4]';
f43=[0, 0, 1]';
f53=[1, 0, 1]';
surface3.contact_pts=[p13 p23 p33 p43 p53];
surface3.contact_fs=[f11 f21 f31 f41 f51];
surface3.mu=0.65;

p14=[4, 22,   2]';
p24=[6, 22,   3]';
p34=[7, 26,   3.5]';
p44=[5, 28,   2.5]';
p54=[4, 26,   2]';
f14=[0, 1, 1]';
f24=[0, 0, 0.8]';
f34=[0.2, 0, 0.4]';
f44=[0, 0, 1]';
f54=[1, 0, 1]';
surface4.contact_pts=[p14 p24 p34 p44 p54];
surface4.contact_fs=[f14 f24 f34 f44 f54];
surface4.mu=0.65;

surfaces=[surface1,surface2,surface3,surface4];
CoM=[0 0 5]';

U_all=zeros(6,4*n_stance);
j=1;
for stance=(1:n_stance)
    if mod(stance, 2)==0
        surf=[surfaces(j),surfaces(j+1)];
        j=j+1;
    else
        surf=(surfaces(j));
    end
    [Wi ,Ui] =  gravito_inertial_wrench(surf, CoM);
    U_all(:,(stance-1)*4+1:stance*4)=Ui;
end


pp_stance = 100;
N = pp_stance;
T = 15;
m = 3.0;
g = [0;0;9.81];
U_stance=zeros(6,4*2);
goal=[ 0, 0,  0  ; 
       4, 7, 14  ;
       5, 5,  5 ];

X=zeros(6,n_stance*N);
sol_all=zeros(9,N);
for j=(1:n_stance-1)
    % we want to shift from a stance to another
    U_stance(:,1:4)=U_all(:,(j-1)*4+1:j*4);
    U_stance(:,5:8)=U_all(:,j*4+1:(j+1)*4);
    
    opti = casadi.Opti();
    
    x   = opti.variable(6,N);
    p   = x(1:3,:);
    dp  = x(4:6,:);
    ddp = opti.variable(3,N); % input
    w   = opti.variable(6,N);
    
    stance1 = opti.parameter(1, N);
    stance2 = opti.parameter(1, N);
    opti.set_value(stance1(1,1),1);
    opti.set_value(stance2(1,1),0);
    
    f = @(dp,ddp) [dp; ddp];
    dt = T/N;
    
    py_des = (0:dt:T);
    
    opti.minimize(trace(ddp'*ddp));
    
    cc = 0;
    
    
    for r_value = (2:pp_stance)
        if r_value < pp_stance/2
           opti.set_value(stance1(1,r_value),1);
           opti.set_value(stance2(1,r_value),0);
        else 
           opti.set_value(stance1(1,r_value),0);
           opti.set_value(stance2(1,r_value),1);
        end
        
        k1 = f(dp(:,r_value),  ddp(:, r_value));
        x_next = x(:,r_value-1) + dt*k1;    
        
        opti.subject_to(x(:,r_value) == x_next);
    
        opti.subject_to(w(1:3, r_value) == m*(g-ddp(:,r_value)));
        opti.subject_to(w(4:6, r_value) == cross(p(:,r_value),m*(g-ddp(:,r_value))));
        
        
        opti.subject_to(dp(:,r_value) >= zeros(3,1))
        
        opti.subject_to(stance1(1,r_value)*U_stance(:,1:4)'*w(:, r_value) <= 0)
        opti.subject_to(stance2(1,r_value)*U_stance(:,5:8)'*w(:, r_value) <= 0)

        
    end
    opti.set_value(stance1(1,N),1);
    opti.set_value(stance2(1,N),0);
    if j==1
        opti.subject_to(x(:,1) == [0;0;5;0;0;0])
        opti.subject_to(x(1:3,N) == [goal(:,j)])
    elseif j==N
        opti.subject_to(x(1:3,1) == [goal(:,j-1)])
        opti.subject_to(x(:,N) == [goal(:,j);0;0;0])
    else
        opti.subject_to(x(1:3,1) == [goal(:,j-1)])
        opti.subject_to(x(1:3,N) == [goal(:,j)])
    end

    opti.solver('ipopt');
    sol = opti.solve();
    
    plot(sol.value(ddp),'o');
    
    sol_all(:,(stance-1)*pp_stance+1:stance*pp_stance)=[sol.value(p);sol.value(dp);sol.value(ddp);];
    X(:,(j-1)*pp_stance+1:j*pp_stance)=[sol.value(p);sol.value(dp)];

end
p_seq = X(1:3,1:N);
plot3(p_seq(1,:),p_seq(2,:),p_seq(3,:))
grid on
ax = gca;
ax.XAxis.Limits = [-10 10];
ax.YAxis.Limits = [0 18];
ax.ZAxis.Limits = [0 10];

