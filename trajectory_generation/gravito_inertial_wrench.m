function [wrench, UGI]= gravito_inertial_wrench(surfaces, CoM)
%   GRAVITO_INERTIAL_WRENCH 
%   Given an array of surfaces, where surface is a structure containing 
%   surface.contact_pts=[point1, point2,...]
%   surface.contact_fs=[force1, force2,...]
%   surface.mu= static friction coefficient of the surface
%   this function plots contact wrenches and computes the gravito-inertial
%   wrench of the stance

    %number of contact planes
    planes=length(surfaces);

    %initialization of useful information
    planes_coef=zeros(planes,4);
    MIN_MAX_x=zeros(2,planes);
    MIN_MAX_y=zeros(2,planes);
    MIN_MAX_z=zeros(2,planes);
    tx_all=zeros(3,planes);
    ty_all=zeros(3,planes);
    t_all=zeros(3,planes);
    tn_all=zeros(3,planes);
    ft_all=zeros(3,planes);
    fn_all=zeros(3,planes);
    Cop_all=zeros(3,planes);
    Fcop_all=zeros(3,planes);
    Tau_cop_all=zeros(3,planes);
    wall = zeros(6*planes,1);
    Astance = zeros(6,6*planes);
    Vall = zeros(6*planes,4*planes);

    %for each contact plane
    for plane=(1:planes)

        %number of contact points
        np=numel(surfaces(plane).contact_pts(1,:));
        
        %find the plane wich contains contact points
        planes_coef(plane,:)=getPlane(surfaces(plane).contact_pts(:,1:3));
        disp(planes_coef(plane,:));

        %plane visualization
        MIN_MAX_x(:,plane)=[min(surfaces(plane).contact_pts(1,:)),max(surfaces(plane).contact_pts(1,:))];
        MIN_MAX_y(:,plane)=[min(surfaces(plane).contact_pts(2,:)),max(surfaces(plane).contact_pts(2,:))];
        MIN_MAX_z(:,plane)=[min(surfaces(plane).contact_pts(3,:)),max(surfaces(plane).contact_pts(3,:))];
        
        [x1 y1] = meshgrid(MIN_MAX_x(1,plane):0.1:MIN_MAX_x(2,plane), MIN_MAX_y(1,plane):0.1:MIN_MAX_y(2,plane)); % Generate x and y data
        z1=(-planes_coef(plane,1)*x1-planes_coef(plane,2)*y1-planes_coef(plane,4))/planes_coef(plane,3); % Solve for z data
        surf(x1,y1,z1) %Plot the surface
        hold on
        
        %contact points visualization 
        px=surfaces(plane).contact_pts(1,:);
        py=surfaces(plane).contact_pts(2,:);
        pz=surfaces(plane).contact_pts(3,:);
        plot3(px,py,pz,'r')
        
        %normal and tangential vectors wrt plane
        tx=[1;0;(-planes_coef(plane,1)-planes_coef(plane,4))/planes_coef(plane,3)];
        ty=[0;1;(-planes_coef(plane,2)-planes_coef(plane,4))/planes_coef(plane,3)];
        t=tx_all(:,plane)+ty_all(:,plane);
        n=planes_coef(plane,1:3)';
        tx_all(:,plane)=tx;
        ty_all(:,plane)=ty;
        t_all(:,plane)=t;
        tn_all(:,plane)=n;
        
        %decomposition of contact forces and creation of wrench
        V=zeros(3*np,4); %span of cone
        X=zeros(np,4);
        Y=zeros(np,4);
        Z=zeros(np,4);
        mu=surfaces(plane).mu;
        Cop=zeros(3,1); %center of pressure 
        Fcop=zeros(3,1);  %Total contact force on the surface
        Tau_cop=zeros(3,1); %Total torque 
        total_sum=0;
        
        %for each contact point
        for i=(1:np)
            fp=surfaces(plane).contact_fs(:,i);
            f_n=n*fp'*n;
            f_tx=tx*fp'*tx;
            f_ty=ty*fp'*ty;
            f_t=f_tx+f_ty;
            fn_all(:,plane)=f_n;
%           Ftx_all(:,plane)=f_tx;
%           Fty_all(:,plane)=f_ty;
            ft_all(:,plane)=f_t;
            f=f_n+f_t;

            %Force visualization 
            quiver3(px(i),py(i),pz(i),fp(1),fp(2),fp(3),'y','LineWidth',2)
            quiver3(px(i),py(i),pz(i),f_n(1),f_n(2),f_n(3),'r','LineWidth',2)
            quiver3(px(i),py(i),pz(i),f_t(1),f_t(2),f_t(3),'g','LineWidth',2)
            quiver3(px(i),py(i),pz(i),f_tx(1),f_tx(2),f_tx(3),'b','LineWidth',2)
            quiver3(px(i),py(i),pz(i),f_ty(1),f_ty(2),f_ty(3),'b','LineWidth',2)
        
            %contact cone visualization
            
            [x,y,z]=cylinder((0:0.1:1)*mu*norm(f_n),4);
        
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
            % Save vectors for span form in surface frame
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
            
            %Save cone base vector
            V(3*(i-1)+1:3*(i-1)+3,1:4)=Vi;
            
            Cop = Cop*total_sum/(total_sum+norm(f)) + norm(f)/(total_sum+norm(f))*[px(i); py(i); pz(i)]; 
            total_sum=total_sum+norm(f);
            
            Fcop= Fcop+fp; 
            tau=cross([px(i), py(i),pz(i)]',-fp);
            Tau_cop=Tau_cop+tau;
            
        end % here we finish to consider all contantact points in a plane
        
        %Save forces, torques and Cop of the surface
        Fcop_all(:,plane)=Fcop;
        Tau_cop_all(:,plane)=Tau_cop;
        Cop_all(:,plane)=Cop;
        
        %face form
        U=zeros(3*np,4);
        for i=(1:3:3*np)
            for j=(1:4)
                if j==4
                    ui=cross(V(i:i+2,1),V(i:i+2,j));
                else
                    ui=cross(V(i:i+2,j+1),V(i:i+2,j));
                end
                % Visualization of face 
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
       
        %Surface wrench creation
        w=[Fcop;Tau_cop];
        Asurf=zeros(6,3*np);
        f_full=zeros(3*np,1);
        j=1;
        for i=(1:3:np*3)
            Asurf(1:3,i:i+2)=eye(3);
            Asurf(4:6,i:i+2)=-skew([px(j),py(j),pz(j)]');
            f_full(i:i+2)=f;
            j=j+1;
        end %end of Asurf
        
        w1=Asurf*f_full;
        
        Vsurf=Asurf*V;
        Astance(:,(plane-1)*6+1:plane*6)=[          -R,   zeros(3,3) ;
                                        -skew(Cop)*R,           -R];


        wall((plane-1)*6+1:plane*6) = w1;
        Vall((plane-1)*6+1:plane*6, (plane-1)*4+1:plane*4) = Vsurf;
    end
    
    % Gravito-Inertial wrench
    wGI=Astance*wall;
    disp(wGI);
    VGI=Astance*Vall;
    VGI_single=zeros(6,4);
    for i=(1:4)
        vec_vgi = zeros(2,3);
        for j = (1:planes)
            vec_vgi(1,:) = vec_vgi(1,:) + [VGI(1,i+(j-1)*4),VGI(2,i+(j-1)*4),VGI(3,i+(j-1)*4)];
            vec_vgi(2,:) = vec_vgi(2,:) + [VGI(4,i+(j-1)*4),VGI(5,i+(j-1)*4),VGI(6,i+(j-1)*4)];
            
        end
        VGI_single(1:3,i) = vec_vgi(1,:)';
        VGI_single(4:6,i) = vec_vgi(2,:)';
        
        quiver3(CoM(1),CoM(2),CoM(3),vec_vgi(1,1),vec_vgi(1,2),vec_vgi(1,3),'g','LineWidth',2)
        
        %quiver3(CoM(1),CoM(2),CoM(3),vec_vgi(2,1),vec_vgi(2,2),vec_vgi(2,3),'r','LineWidth',2)
        
    end
    disp(VGI_single)
    axis equal
    disp(rank(VGI_single))
    wrench=wGI;
    
    UGI=zeros(6,4);
    for j=(1:4)
        if j==4
            ui_f=cross(VGI_single(1:3,1),VGI_single(1:3,j));
            ui_t=cross(VGI_single(4:6,1),VGI_single(4:6,j));
        else
            ui_f=cross(VGI_single(1:3,j+1),VGI_single(1:3,j));
            ui_t=cross(VGI_single(4:6,j+1),VGI_single(4:6,j));
        end
        UGI(1:3,j) = ui_f;
        UGI(4:6,j) = ui_t;
    end %end of U
    UGI
    
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