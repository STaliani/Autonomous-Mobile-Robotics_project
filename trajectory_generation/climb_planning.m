clc
clear all
close all
addpath('casadi-windows-matlabR2016a-v3.5.5')
import casadi.*

%% Constant definition
global ss_time ds_time 

n_stance = 16;                   % number of different conntact configurations
n_surfs = count_surfs(n_stance); % number of surfaces to be generated to have desired stances
t_stance = 1;                    % 0.5 - 2 s per passo 
ss_time= t_stance*0.7;           % time in single support  
ds_time= t_stance*0.3;           % time in double support
T = t_stance*n_stance/2;         % Total simulation time 
dT= 0.1;                         % Integration time
m = 50;                          % 50 kg
g = [0;0;-9.81]; 
eps= 1e-3; 
MIN_X=-34;
MAX_X= 34;

goal= [0; 0.3; 0.40];            % desired final position
%% Get contact surfaces

surfaces= get_surf('hands_tilt',n_surfs);
CoM=[0 0 0.1]'; 

%% Dynamics
X  = casadi.MX.sym('X',6);
U  = casadi.MX.sym('U',3);
p  = X(1:3);
dp = X(4:6);
w  = [m*(g-U);
       cross(p,m*(g-U))]; 
W  = casadi.Function('GI_wrench',{X,U},{w});

X_dot = [dp; U];
F = casadi.Function('continuous_dynamics', {X,U},{X_dot});

% discretization
X_k1 = casadi.MX.sym('X_k1', 6);
U_k1 = casadi.MX.sym('U_k1', 3);

update = X_k1 + dT *  F(X_k1, U_k1); %explicit euler
F_k = casadi.Function('discrete_dynamics', {X_k1, U_k1}, {update});


%% Optimization problem

N = round(T/dT);
opti = casadi.Opti();

Xo = opti.variable(6, N+1);   % state
Uo = opti.variable(3, N+1);   % input

opti.minimize(trace(Uo'*Uo)); % cost function
j=1;
next_var=1;

z = 0.7;

for n = 1:N
   %select stance
   var = find_stance(n*dT);
   if next_var == var(1)
      if mod(var(1), 2)~=0
          surf=[surfaces(j),surfaces(j+1),surfaces(j+2),surfaces(j+3)];
          z = min([mean(surfaces(j).contact_pts(3,:)),mean(surfaces(j+1).contact_pts(3,:)),mean(surfaces(j+2).contact_pts(3,:)),mean(surfaces(j+3).contact_pts(3,:))]) + 0.7;
          j=j+1;
      else
          surf=[surfaces(j),surfaces(j+1),surfaces(j+3)];
          z = min([mean(surfaces(j).contact_pts(3,:)),mean(surfaces(j+1).contact_pts(3,:)),mean(surfaces(j+2).contact_pts(3,:))]) + 0.7;
      end
      Ui =  gravito_inertial_wrench(surf, CoM);
      next_var=next_var+1;
   end
   opti.subject_to(Ui*W(Xo(:,n),Uo(:,n)) <= 0);        % constrain the giw 
   opti.subject_to(Xo(:,n+1) == F_k(Xo(:,n),Uo(:,n))); % impose system dynamics
  
   opti.subject_to(Xo(3,n) >= z-0.1);                  % change 0.5 on the base
                                                       % of the min elevation
   opti.subject_to(Xo(3,n) <= z);                      % impose almost constant CoM height
end

%%impose initial position
opti.subject_to(Xo(1,1) == 0)
opti.subject_to(Xo(2,1) == 0)
opti.subject_to(Xo(3,1) == CoM(3))
opti.subject_to(Xo(4,1) == 0)
opti.subject_to(Xo(5,1) == 0)
opti.subject_to(Xo(6,1) == 0)
%%impose final position
opti.subject_to(Xo(2,N) == goal(2))
opti.subject_to(Xo(1,N) == goal(1))
opti.subject_to(Xo(4,N) == 0)
opti.subject_to(Xo(5,N) == 0)
opti.subject_to(Xo(6,N) == 0)

%solve problem
opti.solver('ipopt');
sol = opti.solve();

% To have the average through 30 attemts uncomment the following for
% time = 0;
% for i = 1:30
%     sol = opti.solve();
%     time = time+sol.stats.t_wall_total/30;
% end

%% solution visualization 
res = sol.value(Xo);
p_seq = res(1:3,1:N);
plt_seq(p_seq,N,t_stance,dT)

%% solution visualization 
res = sol.value(Xo);
inputs = sol.value(Uo);
p_seq = res(1:3,1:N);
plt_seq(p_seq,N,t_stance,dT)
axis equal

figure Name 'Cartesian Position'
plot(0:dT:T, res(1,:), 'LineWidth', 2, 'Color', 'c')
hold on
plot(0:dT:T, res(2,:), 'LineWidth', 2,'LineStyle','--')
plot(0:dT:T, res(3,:), 'LineWidth', 2)
grid on
xlabel('time [s]');
ylabel('position [m]');
legend('$X$','$Y$','$Z$', 'interpreter','latex','FontSize',14)

figure Name 'Cartesian Velocity'
plot(0:dT:T, res(4,:), 'LineWidth', 2, 'Color', 'c')
hold on
plot(0:dT:T, res(5,:), 'LineWidth', 2,'LineStyle','--')
plot(0:dT:T, res(6,:), 'LineWidth', 2)
grid on
xlabel('time [s]');
ylabel('velocity [m/s]');
legend('$\dot{X}$','$\dot{Y}$','$\dot{Z}$', 'interpreter','latex','FontSize',14)

figure Name 'Commanded Accelerations'
plot(0:dT:T, inputs(1,:), 'LineWidth', 2, 'Color', 'c')
hold on
plot(0:dT:T, inputs(2,:), 'LineWidth', 2,'LineStyle','--')
plot(0:dT:T, inputs(3,:), 'LineWidth', 2)
grid on
xlabel('time [s]');
ylabel('acceleration [m/s^2]');
legend('$\ddot{X}$','$\ddot{Y}$','$\ddot{Z}$', 'interpreter','latex','FontSize',14)



%% Auxiliary functions

function stance = find_stance(t) 
    % n_support defines the number of supports we are using 
    % n_stance defines in which stance we are
    global ss_time ds_time
    n_support=1;
    res = mod(t,ss_time+ds_time);
    if res>ss_time
        n_support=2;
    end
    n_stance = fix(t/(ss_time+ds_time));
    stance = [n_stance*2+n_support, n_support]; 
end

function surfs = count_surfs(n_stance)
    if or(n_stance==1,n_stance==2)
        surfs=4;
    else
        surfs= round((n_stance+7)/2);
    end
end

function a = plt_seq(p_seq,N,t_stance,dT)
    global ss_time ds_time
    for i= 1:N
        s = mod(i-1,t_stance/dT);
        if s<=ds_time/dT
                plot3(p_seq(1,i),p_seq(2,i),p_seq(3,i),'o','Color','r')

                hold on;
        else
                plot3(p_seq(1,i),p_seq(2,i),p_seq(3,i),'o','Color','g');
        
        end
    end
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    a=1;
end

