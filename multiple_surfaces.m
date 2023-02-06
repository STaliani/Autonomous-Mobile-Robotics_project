clc
clear all 

%% robot characteristics
mass=0.1; %kg
g=[0,0,-9.81]; %gravity
fp=mass*g;  
CoM=[0,0,10]'; %position of CoM

%points of contact
p1=[ 4,  2,  2;
    -4, -8,  4];
p2=[ 6,  2,  3;
    -6, -8,  4];
p3=[ 7,  6,  3.5;
    -7, -4,  2];
p4=[ 5,  8, 2.5;
    -5, -2,  1];
p5=[ 4,  6,  2;
    -4, -4,  2];
px=[p1(1,1),p2(1,1),p3(1,1),p4(1,1),p5(1,1),p1(1,1);
    p1(2,1),p2(2,1),p3(2,1),p4(2,1),p5(2,1),p1(2,1)];
py=[p1(1,2),p2(1,2),p3(1,2),p4(1,2),p5(1,2),p1(1,2);
    p1(2,2),p2(2,2),p3(2,2),p4(2,2),p5(2,2),p1(2,2)];
pz=[p1(1,3),p2(1,3),p3(1,3),p4(1,3),p5(1,3),p1(1,3);
    p1(2,3),p2(2,3),p3(2,3),p4(2,3),p5(2,3),p1(2,3)];
p1=p1';
p2=p2';
p3=p3';
p4=p4';
p5=p5';
np=numel(px)-1;
%% define contact surface
[x1 y1] = meshgrid(0:0.1:10); % Generate x and y data
z1 = 0.5*x1 + 0*y1 ; % Solve for z data
surf(x1,y1,z1) %Plot the surface
hold on
[x2 y2] = meshgrid(-10:0.1:0); % Generate x and y data
z2 = 0*x2 - 0.5*y2 ; % Solve for z data
surf(x2,y2,z2) %Plot the surface
ax = gca;
axis equal
mu=0.65; % coeff, gomma su cemento asciutto

plot3(px(1,:),py(1,:),pz(1,:),"r");
plot3(px(2,:),py(2,:),pz(2,:),"r");

planes=2;
alpha_tan = [0, 0.5]; % rotation w.r.t. x
beta_tan = [-0.5, 0]; % rotation w.r.t. y

%directions 
tx=[ 1, 0, 0.5;
    -1, 0,   0]';

ty=[0,  1,   0;
    0, -1, 0.5]';

n=zeros(3,planes);
t=zeros(3,planes);
for dir=(1:2)
    tx(:,dir)=tx(:,dir)/norm(tx(:,dir));
    ty(:,dir)=ty(:,dir)/norm(ty(:,dir));
    t(:,dir)=tx(:,dir)+ty(:,dir);
    t(:,dir)=t(:,dir)/norm(t(:,dir));
    n(:,dir)=cross(tx(:,dir), ty(:,dir));
    n(:,dir)=n(:,dir)/norm(n(:,dir));
end
%% define contact forces
Fn_all=zeros(3,planes);
Ftx_all=zeros(3,planes);
Fty_all=zeros(3,planes);
Ft_all=zeros(3,planes);
Cop_all=zeros(3,planes);

wall = zeros(6*planes,1);
Astance = zeros(6,6*planes);
Vall = zeros(6*planes,4*planes);

for plane=(1:planes)
    
    f_n=n(:,plane)*fp*n(:,plane);
    f_tx=tx(:,plane)*fp*tx(:,plane);
    f_ty=ty(:,plane)*fp*ty(:,plane);
    f_t=f_tx+f_ty;

    Fn_all(:,plane)=f_n;
    Ftx_all(:,plane)=f_tx;
    Fty_all(:,plane)=f_ty;
    Ft_all(:,plane)=f_t;
    
    if norm(f_t)<=mu*f_n 
        static=true;
    else
        static=false;
    end
    f=f_n+f_t;

    V=zeros(3*np,4); %span of cone
    X=zeros(np,4);
    Y=zeros(np,4);
    Z=zeros(np,4);
    for i=(1:5)
        quiver3(px(plane,i),py(plane,i),pz(plane,i),f(1),f(2),f(3),'m','LineWidth',2)
        quiver3(px(plane,i),py(plane,i),pz(plane,i),fp(1),fp(2),fp(3),'y','LineWidth',2)
        quiver3(px(plane,i),py(plane,i),pz(plane,i),-f_n(1),-f_n(2),-f_n(3),'r','LineWidth',2)
        quiver3(px(plane,i),py(plane,i),pz(plane,i),f_t(1),f_t(2),f_t(3),'g','LineWidth',2)
        quiver3(px(plane,i),py(plane,i),pz(plane,i),f_tx(1),f_tx(2),f_tx(3),'b','LineWidth',2)
        quiver3(px(plane,i),py(plane,i),pz(plane,i),f_ty(1),f_ty(2),f_ty(3),'b','LineWidth',2)
    
        %wrench cone
        [x,y,z]=cylinder((0:0.1:1)*mu*norm(f_n),4);
    
        % Rotate the cone to align with the direction vector
        alpha=atan(alpha_tan(plane));
        beta = atan(beta_tan(plane));
        Rx = [        1,           0,                0 ;
                      0,  cos(alpha),       sin(alpha) ;
                      0, -sin(alpha),       cos(alpha)];
        Ry = [cos(beta),        0,        sin(beta) ;
                      0,        1,                0 ;
             -sin(beta),        0,        cos(beta)];   
        R = Ry*Rx;
        
        for j = 1:size(x, 2)
            [x(:, j), y(:, j), z(:, j)] = apply_rotation(x(:, j), y(:, j), z(:, j), R);
        end
    
        % Translate the cone to the origin
        x = x + px(plane,i);
        y = y + py(plane,i);
        z = z + pz(plane,i);
        temp1=R'*[x(1,1:4);y(1,1:4);z(1,1:4)];
        temp =R'*[x(11,1:4);y(11,1:4);z(11,1:4)];
        X(i,1:4)=temp(1,1:4);
        Y(i,1:4)=temp(2,1:4);
        Z(i,1:4)=temp(3,1:4);
        % Plot the cone
        surf(x, y, z);
        Vi=zeros(3,4);
        for j=(1:4)
            v=[X(i,j)-temp1(1,j),Y(i,j)-temp1(2,j),Z(i,j)-temp1(3,j)]';
            Vi(1:3,j)=v;
        end
        %[Vi,R, X(i,:),Y(i,:),Z(i,:)]=create_cone(mu, f_n, alpha_tan(plane),beta_tan(plane),[px(plane,i),py(plane,i),pz(plane,i)]);
        V(3*(i-1)+1:3*(i-1)+3,1:4)=Vi;
    end % here we finish to consider all contantact points in a plane

    U=zeros(3*np,4);
    for i=(1:3:3*np)
        for j=(1:4)
            if j==4
                ui=cross(V(i:i+2,1),V(i:i+2,j));
            else
                ui=cross(V(i:i+2,j+1),V(i:i+2,j));
            end
            temp1=R*[X((i+2)/3,j);Y((i+2)/3,j);Z((i+2)/3,j)];
            temp3=R*ui;
            if j==4            
                temp2=R*[X((i+2)/3,1);Y((i+2)/3,1);Z((i+2)/3,1)];
                quiver3((temp1(1)+temp2(1))/2,(temp1(2)+temp2(2))/2,(temp1(3)+temp2(3))/2,temp3(1),temp3(2),temp3(3),'b','LineWidth',2)
            else
                temp2=R*[X((i+2)/3,j+1);Y((i+2)/3,j+1);Z((i+2)/3,j+1)];
                quiver3((temp1(1)+temp2(1))/2,(temp1(2)+temp2(2))/2,(temp1(3)+temp2(3))/2,temp3(1),temp3(2),temp3(3),'b','LineWidth',2)
            end
            U(i:i+2,j)=ui;
        end
    end %end of U
    
    
    
    Cop = [mean(px(plane,1:5)), mean(py(plane,1:5)), mean(pz(plane,1:5))];
    Cop_all(:,plane)=Cop';
    Fcop= f*5;
    tau=[0 0 0]';
    Tau_cop=[0 0 0]';
    for i=(1:5)
        tau=cross([px(i), py(i),pz(i)]',-f);
        Tau_cop=Tau_cop+tau;
    end %end Tau_cop
    
    w=[Fcop;Tau_cop];
    Asurf=zeros(6,3*np);
    f_full=zeros(3*np,1);
    j=1;
    for i=(1:3:15)
        Asurf(1:3,i:i+2)=eye(3);
        Asurf(4:6,i:i+2)=-skew([px(j),py(j),pz(j)]');
        f_full(i:i+2)=f;
        j=j+1;
    end %end of Asurf
    
    w1=Asurf*f_full
    
    Vsurf=Asurf*V;
       
    Astance(:,(plane-1)*6+1:plane*6)=[          -R,   zeros(3,3) ;
                                        -skew(Cop)*R,           -R];
    wall((plane-1)*6+1:plane*6) = w1;
    Vall((plane-1)*6+1:plane*6, (plane-1)*4+1:plane*4) = Vsurf;

end %end of planes


wGI=Astance*wall

VGI=Astance*Vall

for i=(1:4)
    vec_vgi = zeros(2,3);
    for j = (1:planes)
        vec_vgi(1,:) = vec_vgi(1,:) + [-VGI(1,i+(j-1)*4),-VGI(2,i+(j-1)*4),-VGI(3,i+(j-1)*4)];
        vec_vgi(2,:) = vec_vgi(2,:) + [-VGI(4,i+(j-1)*4),-VGI(5,i+(j-1)*4),-VGI(6,i+(j-1)*4)];
    end
    quiver3(CoM(1),CoM(2),CoM(3),vec_vgi(1,1),vec_vgi(1,2),vec_vgi(1,3),'g','LineWidth',2)
    
    %quiver3(CoM(1),CoM(2),CoM(3),vec_vgi(2,1),vec_vgi(2,2),vec_vgi(2,3),'r','LineWidth',2)
    
end
%create_cone(mu,vec_vgi(1,:),atan2(vec_vgi(1,3),vec_vgi(1,1)),atan2(vec_vgi(1,3),vec_vgi(1,2)),CoM)
%create_cone(1,vec_vgi(2,:),atan2(vec_vgi(2,3),vec_vgi(2,1)),atan2(vec_vgi(2,3),vec_vgi(2,2)),CoM)
axis equal







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
function [Vi,R,X,Y,Z] = create_cone(mu, f_n, alpha_tan,beta_tan,p)
        %wrench cone
        [x,y,z]=cylinder((0:0.1:1)*mu*norm(f_n),4);
    
        % Rotate the cone to align with the direction vector
        alpha=atan(alpha_tan);
        beta = atan(beta_tan);
        Rx = [        1,           0,                0 ;
                      0,  cos(alpha),       sin(alpha) ;
                      0, -sin(alpha),       cos(alpha)];
        Ry = [cos(beta),        0,        sin(beta) ;
                      0,        1,                0 ;
             -sin(beta),        0,        cos(beta)];   
        R = Ry*Rx;
        
        for j = 1:size(x, 2)
            [x(:, j), y(:, j), z(:, j)] = apply_rotation(x(:, j), y(:, j), z(:, j), R);
        end
    
        % Translate the cone to the origin
        x = x + p(1);
        y = y + p(2);
        z = z + p(3);
        temp1=R'*[x(1,1:4);y(1,1:4);z(1,1:4)];
        temp =R'*[x(11,1:4);y(11,1:4);z(11,1:4)];
        X(i,1:4)=temp(1,1:4);
        Y(i,1:4)=temp(2,1:4);
        Z(i,1:4)=temp(3,1:4);
        % Plot the cone
        surf(x, y, z);
        Vi=zeros(3,4);
        for j=(1:4)
            v=[X(i,j)-temp1(1,j),Y(i,j)-temp1(2,j),Z(i,j)-temp1(3,j)]';
            Vi(1:3,j)=v;
        end

end

