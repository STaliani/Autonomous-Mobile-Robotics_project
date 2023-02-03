clc
clear all 

%% robot characteristics
mass=0.1; %kg
g=[0,0,-9.81]; %gravity
fp=mass*g;

%points of contact
p1=[4, 2,   2]';
p2=[6, 2,   3]';
p3=[6, 6,   3]';
p4=[5, 8, 2.5]';
p5=[3, 6, 1.5]';
px=[p1(1),p2(1),p3(1),p4(1),p5(1),p1(1)];
py=[p1(2),p2(2),p3(2),p4(2),p5(2),p1(2)];
pz=[p1(3),p2(3),p3(3),p4(3),p5(3),p1(3)];
np=numel(px)-1;
%% define contact surface
[x y] = meshgrid(0:0.1:10); % Generate x and y data
z = 0.5*x + 0*y ; % Solve for z data
surf(x,y,z) %Plot the surface
hold on
ax = gca;
ax.XAxis.Limits = [0 10];
ax.YAxis.Limits = [0 10];
ax.ZAxis.Limits = [0 10];
mu=0.65; % coeff, gomma su cemento asciutto

plot3(px,py,pz,"r");


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
%tau=cross(p,f); %Nm contact torque

%w=[f; tau];
V=zeros(3*np,4);
X=zeros(np,4);
Y=zeros(np,4);
Z=zeros(np,4);
for i=(1:5)
    quiver3(px(i),py(i),pz(i),f(1),f(2),f(3),'m','LineWidth',2)
    quiver3(px(i),py(i),pz(i),fp(1),fp(2),fp(3),'y','LineWidth',2)
    quiver3(px(i),py(i),pz(i),-f_n(1),-f_n(2),-f_n(3),'r','LineWidth',2)
    quiver3(px(i),py(i),pz(i),f_t(1),f_t(2),f_t(3),'g','LineWidth',2)
    quiver3(px(i),py(i),pz(i),f_tx(1),f_tx(2),f_tx(3),'b','LineWidth',2)
    quiver3(px(i),py(i),pz(i),f_ty(1),f_ty(2),f_ty(3),'b','LineWidth',2)

    %wrench cone
    [x,y,z]=cylinder((0:0.1:1)*mu*norm(f_n),4);

    % Rotate the cone to align with the direction vector
    alpha=atan(-0.5);
    R = [cos(alpha),        0,       sin(alpha) ;
                  0,        1,                0 ;
        -sin(alpha),        0,       cos(alpha)];                                                  
    for j = 1:size(x, 2)
        [x(:, j), y(:, j), z(:, j)] = apply_rotation(x(:, j), y(:, j), z(:, j), R);
    end

    % Translate the cone to the origin
    x = x + px(i);
    y = y + py(i);
    z = z + pz(i);
    X(i,1:4)=x(11,1:4);
    Y(i,1:4)=y(11,1:4);
    Z(i,1:4)=z(11,1:4);
    % Plot the cone
    h = surf(x, y, z);




    Vi=zeros(3,4);
    for j=(1:4)
        v=[X(i,j)-x(1,j),Y(i,j)-y(1,j),Z(i,j)-z(1,j)]';
        Vi(1:3,j)=v;
    end
    V(3*(i-1)+1:3*(i-1)+3,1:4)=Vi;
end

U=zeros(3*np,4);
for i=(1:3:3*np)
    for j=(1:4)
        if j==4
            ui=cross(V(i:i+2,1),V(i:i+2,j));
        else
            ui=cross(V(i:i+2,j+1),V(i:i+2,j));
        end
        if j==4
            quiver3((X((i+2)/3,j)+X((i+2)/3,1))/2,(Y((i+2)/3,j)+Y((i+2)/3,1))/2,(Z((i+2)/3,j)+Z((i+2)/3,1))/2,ui(1),ui(2),ui(3),'b','LineWidth',2)
        else
            quiver3((X((i+2)/3,j)+X((i+2)/3,j+1))/2,(Y((i+2)/3,j)+Y((i+2)/3,j+1))/2,(Z((i+2)/3,j)+Z((i+2)/3,j+1))/2,ui(1),ui(2),ui(3),'b','LineWidth',2)
        end
        U(i:i+2,j)=ui;
    end
end



Cop= [mean(px(1:5)), mean(py(1:5)), mean(pz(1:5))];
Fcop= f*5;
tau=[0 0 0]';
Tau_cop=[0 0 0]';
for i=(1:5)
    tau=cross([px(i), py(i),pz(i)]',-f);
    Tau_cop=Tau_cop+tau;
end

w=[Fcop;Tau_cop];
Asurf=zeros(6,3*np);
f_full=zeros(3*np,1);
j=1;
for i=(1:3:15)
    Asurf(1:3,i:i+2)=eye(3);
    Asurf(4:6,i:i+2)=-skew([px(j),py(j),pz(j)]');
    f_full(i:i+2)=f;
    j=j+1;
end

w1=Asurf*f_full

Vsurf=Asurf*V;


Astance=[    -R,       zeros(3,3) ;
         -skew(Cop)*R,         -R];

wGI=Astance*w1

VGI=Astance*Vsurf

for i=(1:4)
    quiver3(Cop(1),Cop(2),Cop(3),-VGI(1,i),-VGI(2,i),-VGI(3,i),'g','LineWidth',2)
end










function s=skew(p)
    s=[   0, -p(3),  p(2);
       p(3),     0, -p(1);
      -p(2),  p(1),    0];
end
        
    % Apply rotation to each point
function [x, y, z] = apply_rotation(x, y, z, R)
    points = [x(:), y(:), z(:)];
    rotated_points = (R * points')';
    x = reshape(rotated_points(:, 1), size(x));
    y = reshape(rotated_points(:, 2), size(y));
    z = reshape(rotated_points(:, 3), size(z));
end
