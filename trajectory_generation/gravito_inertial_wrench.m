function UGI = gravito_inertial_wrench(surfaces, CoM)
%   GRAVITO_INERTIAL_WRENCH 
%   Given an array of surfaces, where surface is a structure containing 
%   surface.contact_pts=[point1, point2,...]
%   surface.mu= static friction coefficient of the surface
%   this function plots contact wrenches and computes the gravito-inertial
%   wrench of the stance
    addpath('polytopes_2017_10_04_v1.9\')
    %number of contact planes
    planes=length(surfaces);

    %initialization of useful information
    planes_coef=zeros(planes,4);
    MIN_MAX_x=zeros(2,planes);
    MIN_MAX_y=zeros(2,planes);
    MIN_MAX_z=zeros(2,planes);
    Astance = zeros(6,6*planes);
    Vall = [];

    %for each contact plane
    for plane=(1:planes)

        %number of contact points
        np=numel(surfaces(plane).contact_pts(1,:));
        
        %find the plane wich contains contact points
        planes_coef(plane,:)=getPlane(surfaces(plane).contact_pts(:,1:3));

        %plane visualization
        MIN_MAX_x(:,plane)=[min(surfaces(plane).contact_pts(1,:))-0.1,max(surfaces(plane).contact_pts(1,:))+0.1];
        MIN_MAX_y(:,plane)=[min(surfaces(plane).contact_pts(2,:))-0.1,max(surfaces(plane).contact_pts(2,:))+0.1];
        MIN_MAX_z(:,plane)=[min(surfaces(plane).contact_pts(3,:))-0.1,max(surfaces(plane).contact_pts(3,:))+0.1];
        
        [x1, y1] = meshgrid(MIN_MAX_x(1,plane):0.1:MIN_MAX_x(2,plane), MIN_MAX_y(1,plane):0.1:MIN_MAX_y(2,plane)); % Generate x and y data
        z1=(-planes_coef(plane,1)*x1-planes_coef(plane,2)*y1-planes_coef(plane,4))/planes_coef(plane,3); % Solve for z data
        surf(x1,y1,z1) %Plot the surface
        hold on
        
        %contact points visualization 
        px=surfaces(plane).contact_pts(1,:);
        py=surfaces(plane).contact_pts(2,:);
        pz=surfaces(plane).contact_pts(3,:);
        plot3(px,py,pz,'r')
        
        %normal and tangential vectors wrt plane
        tx=[1;0;(-planes_coef(plane,1))/planes_coef(plane,3)];
        ty=[0;1;(-planes_coef(plane,2))/planes_coef(plane,3)];
        t=tx+ty;
        n=planes_coef(plane,1:3)';
          
        %friction coefficient of the contact suface for cone visualization
        mu=surfaces(plane).mu;

        %for each contact point
        for i=(1:np)
            
            %plane's vector visualization 
             %quiver3(px(i),py(i),pz(i),n(1),n(2),n(3),'r','LineWidth',2)
             %quiver3(px(i),py(i),pz(i),t(1),t(2),t(3),'g','LineWidth',2)
%             quiver3(px(i),py(i),pz(i),tx(1),tx(2),tx(3),'b','LineWidth',2)
%             quiver3(px(i),py(i),pz(i),ty(1),ty(2),ty(3),'b','LineWidth',2)
        
            %friction cone's point generation
            [x,y,z]=cylinder((0:0.1:1)*mu,6);
            x=x*0.1;
            y=y*0.1;
            z=z*0.1;
            % Rotate the cone to align with the direction vector
            alpha= atan(planes_coef(plane,2)/planes_coef(plane,3));
            beta = atan(planes_coef(plane,1)/planes_coef(plane,3));
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
            x = x + px(i);
            y = y + py(i);
            z = z + pz(i);
            
            % Plot the cone
            surf(x, y, z); 
        end % here we finish to consider all contantact points in a plane
        
        % generation of wrench cone's span form of the surface 
        Vsurf=face_span_gen(surfaces(plane).X, surfaces(plane).Y, mu);
            
        %Surface wrench creation
        Cop=Cop_function(surfaces(plane).contact_pts); %center of pressure 
        

        Astance(:,(plane-1)*6+1:plane*6)=[          -R,   zeros(3,3) ;
                                          -skew(Cop)*R,           -R];

        [row1, col1] = size(Vall); 
        [row2, col2] = size(Vsurf);
        Vall = [             Vall, zeros(row1,col2) ;
                 zeros(row2,col1),           Vsurf ];
    end % here we finish operations on single surfaces
    
    % Gravito-Inertial wrench cone
    % span form
    VGI=Astance*Vall;
    s=size(VGI);

    % GIWC span vectors visualization
    for i= 1 : s(2) 
        %quiver3(CoM(1),CoM(2),CoM(3),VGI(1,i),VGI(2,i),VGI(3,i),'g','LineWidth',2)    
        %quiver3(CoM(1),CoM(2),CoM(3),VGI(4,i)/10,VGI(5,i)/10,VGI(6,i)/10,'r','LineWidth',2)
    end

    %VGI=licols(VGI)
    UGI = face_of_span(VGI);

end




%% Complementary functios

function varout = getPlane(xyz)
    x1 = xyz(1,1);
    y1 = xyz(2,1);
    z1 = xyz(3,1);
    x2 = xyz(1,2);
    y2 = xyz(2,2);
    z2 = xyz(3,2);
    x3 = xyz(1,3);
    y3 = xyz(2,3);
    z3 = xyz(3,3);
    
    % Calculate the normal vector of the plane
    v1 = [x2-x1, y2-y1, z2-z1];
    v2 = [x3-x1, y3-y1, z3-z1];
    n = cross(v1, v2);
    
    % Normalize the normal vector
    n = n / norm(n);
    
    % Get the coefficients a, b, c, and d
    a = n(1);
    b = n(2);
    c = n(3);
    d = -(a*x1 + b*y1 + c*z1);
    varout = [a,b,c,d];
    if c<0
        varout=-varout;
    end
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

% function that computes the center of pressure in a surface
function Cop = Cop_function(points)
    np=numel(points(1));
    x=0;
    y=0;
    z=0;
    for i = 1:np
        x=x + points(1,i)/np;
        y=y + points(2,i)/np;
        z=z + points(3,i)/np;
    end
    Cop = [x;y;z];
end
function [Xsub,idx]=licols(X,tol)
%Extract a linearly independent set of columns of a given matrix X
%
%    [Xsub,idx]=licols(X)
%
%in:
%
%  X: The given input matrix
%  tol: A rank estimation tolerance. Default=1e-10
%
%out:
%
% Xsub: The extracted columns of X
% idx:  The indices (into X) of the extracted columns
   if ~nnz(X) %X has no non-zeros and hence no independent columns
       
       Xsub=[]; idx=[];
       return
   end
   if nargin<2, tol=1e-10; end
   
           
     [Q, R, E] = qr(X,0); 
     
     if ~isvector(R)
      diagr = abs(diag(R));
     else
      diagr = abs(R(1));   
     end
     %Rank estimation
     r = find(diagr >= tol*diagr(1), 1, 'last'); %rank estimation
     idx=sort(E(1:r));
     Xsub=X(:,idx);
end