clc
clear all
addpath('casadi-windows-matlabR2016a-v3.5.5')
%addpath('C:/Users/anto1/Documents/MATLAB/casadi-windows-matlabR2016a-v3.5.5')
%addpath('C:/Users/saver/Documents/MATLAB/casadi-windows-matlabR2016a-v3.5.5')
import casadi.*

n_stance = 4;

p11=[-4, -8,   0]';
p21=[-6, -8,   0]';
p31=[-7, -4,   0]';
p41=[-5, -2,   0]';
p51=[-4, -4,   0]';
f11=[0, 1, 1]';
f21=[0, 0, 0.8]';
f31=[0.2, 0, 0.4]';
f41=[0, 0, 1]';
f51=[1, 0, 1]';
surface1.contact_pts=[p11 p21 p31 p41 p51];
surface1.contact_fs=[f11 f21 f31 f41 f51];
surface1.mu=0.65;
p12=[4, 2,  0]';
p22=[6, 2,  0]';
p32=[7, 6,  0]';
p42=[5, 8,  0]';
p52=[4, 6,  0]';
f12=[0, 0, 1]';
f22=[0, 0.1, 0.8]';
f32=[0, 0.2, 0.4]';
f42=[0, 0, 1]';
f52=[0.3, 0, 1]';
surface2.contact_pts=[p12 p22 p32 p42 p52];
surface2.contact_fs=[f12 f22 f32 f42 f52];
surface2.mu=0.65;


p13=[-4, 12,   0]';
p23=[-6, 12,   0]';
p33=[-7, 16,   0]';
p43=[-5, 18,   0]';
p53=[-4, 16,   0]';
f13=[0, 1, 1]';
f23=[0, 0, 0.8]';
f33=[0.2, 0, 0.4]';
f43=[0, 0, 1]';
f53=[1, 0, 1]';
surface3.contact_pts=[p13 p23 p33 p43 p53];
surface3.contact_fs=[f11 f21 f31 f41 f51];
surface3.mu=0.65;

p14=[4, 22,   0]';
p24=[6, 22,   0]';
p34=[7, 26,   0]';
p44=[5, 28,   0]';
p54=[4, 26,   0]';
f14=[0, 1, 1]';
f24=[0, 0, 0.8]';
f34=[0.2, 0, 0.4]';
f44=[0, 0, 1]';
f54=[1, 0, 1]';
surface4.contact_pts=[p14 p24 p34 p44 p54];
surface4.contact_fs=[f14 f24 f34 f44 f54];
surface4.mu=0.65;

surfaces=[surface1,surface2,surface3,surface4];
CoM=[0 0 5]'

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
N = n_stance*pp_stance;
T = 15;
m = 3.0;
g = [0;0;9.81];
opti = casadi.Opti();

x   = opti.variable(6,N+1);
p   = x(1:3,:);
dp  = x(4:6,:);
ddp = opti.variable(3,N+1); % input
w   = opti.variable(6,N+1);

stances_all = opti.parameter(n_stance, N+1);

f = @(dp,ddp) [dp; ddp];
dt = T/N;

py_des = (0:dt:T);

opti.minimize(trace(ddp'*ddp));

for h = (1:n_stance)
   if h == 1
       opti.set_value(stances_all(h,1),1);
   else
       opti.set_value(stances_all(h,1),0);
   end
end

cc = 0;
cont=0;
for r_value = (2:n_stance*pp_stance)
       var=int16(floor(r_value/pp_stance)); 
       
       for h = (1:n_stance)
           if h == var
               opti.set_value(stances_all(h,r_value),1);
           else
               opti.set_value(stances_all(h,r_value),0);
           end
       end
    
    k1 = f(dp(:,r_value),  ddp(:, r_value));
    x_next = x(:,r_value-1) + dt*k1;    
    
    opti.subject_to(x(:,r_value) == x_next);

    opti.subject_to(w(1:3, r_value) == m*(g-ddp(:,r_value)));
    opti.subject_to(w(4:6, r_value) == cross(p(:,r_value),m*(g-ddp(:,r_value))));
    
    %if cc == 30 
    %    opti.subject_to(p == [0;py_des(r_value);5])
    %    cc = 0;
    %end
    %cc = cc +1;
    
    opti.subject_to(dp(:,r_value) >= zeros(3,1))
    
    
    for h = (1:n_stance)
        Ui = U_all(:,(h-1)*4+1:h*4);
        opti.subject_to(stances_all(h,r_value)*Ui'*w(:, r_value) <= 0)
    end
    
end

for h = (1:n_stance)
   if h == n_stance
       opti.set_value(stances_all(h,N+1),1);
   else
       opti.set_value(stances_all(h,N+1),0);
   end
end

opti.subject_to(x(:,1) == [0;0;5;0;0;0])
opti.subject_to(x(1:3,N) == [0;14;5])
opti.subject_to(x(:,N+1) == [0;14;5;0;0;0])

opti.solver('ipopt');
sol = opti.solve();

% plot(sol.value(ddp),'o');

res = sol.value(p);
p_seq = res(:,1:N);

% plot((0:dt:T-dt), p_seq(1,:))
% plot((0:dt:T-dt), p_seq(2,:))
% plot((0:dt:T-dt), p_seq(3,:))

plot3(p_seq(1,:),p_seq(2,:),p_seq(3,:))
grid on
ax = gca;
ax.XAxis.Limits = [-10 10];
ax.YAxis.Limits = [0 18];
ax.ZAxis.Limits = [0 10];



