clc
clear all 

%% robot characteristics
mass=0.1; %kg
g=[0,0,-9.81]; %gravity
fp=mass*g;

%point of contact
p=[1,1, 0.5]';

%% define contact surface
[x y] = meshgrid(0:0.1:2); % Generate x and y data
z = 0.5*x + 0*y ; % Solve for z data
surf(x,y,z) %Plot the surface
hold on
ax = gca;
ax.XAxis.Limits = [-1 2];
ax.YAxis.Limits = [-1 2];
ax.ZAxis.Limits = [-1 2];
mu=0.65; % coeff, gomma su cemento asciutto

%% define contact forces

%directions 
tx=[1, 0,0.5]';
tx=tx/norm(tx);
ty=[0, 1,  0]';
ty=ty/norm(ty);
t=tx+ty;
t=t/norm(t);
n=cross(tx, ty);
n=n/norm(n);

f_n=n*fp*n;
f_tx=tx*fp*tx;
f_ty=ty*fp*ty;
f_t=f_tx+f_ty;

if norm(f_t)<=mu*f_n 
    static=true;
else
    static=false;
end
f=f_n+f_t;
tau=cross(p,f); %Nm contact torque

%w=[f; tau];
quiver3(p(1),p(2),p(3),f(1),f(2),f(3),'m','LineWidth',2)
quiver3(p(1),p(2),p(3),fp(1),fp(2),fp(3),'y','LineWidth',2)
quiver3(p(1),p(2),p(3),-f_n(1),-f_n(2),-f_n(3),'r','LineWidth',2)
quiver3(p(1),p(2),p(3),f_t(1),f_t(2),f_t(3),'g','LineWidth',2)
quiver3(p(1),p(2),p(3),f_tx(1),f_tx(2),f_tx(3),'b','LineWidth',2)
quiver3(p(1),p(2),p(3),f_ty(1),f_ty(2),f_ty(3),'b','LineWidth',2)
quiver3(p(1),p(2),p(3),tau(1),tau(2),tau(3),'g','LineWidth',2)
%plot(w)
%wrench cone
[x,y,z]=cylinder((0:0.1:1)*mu*norm(f_n),4);

% Rotate the cone to align with the direction vector
alpha=atan(-0.5);
R = [cos(alpha),        0,       sin(alpha) ;
              0,        1,                0 ;
    -sin(alpha),        0,       cos(alpha)];                                                  
for i = 1:size(x, 2)
    [x(:, i), y(:, i), z(:, i)] = apply_rotation(x(:, i), y(:, i), z(:, i), R);
end

% Translate the cone to the origin
x = x + p(1);
y = y + p(2);
z = z + p(3);

% Plot the cone
h = surf(x, y, z);

% Apply rotation to each point
function [x, y, z] = apply_rotation(x, y, z, R)
    points = [x(:), y(:), z(:)];
    rotated_points = (R * points')';
    x = reshape(rotated_points(:, 1), size(x));
    y = reshape(rotated_points(:, 2), size(y));
    z = reshape(rotated_points(:, 3), size(z));
end
