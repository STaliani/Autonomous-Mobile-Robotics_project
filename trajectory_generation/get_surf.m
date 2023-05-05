% Select the type of surface: flat, inclined
function surfs = get_surf(type_of_surface, n_surfs)

    switch type_of_surface
        case 'stairs'
            surfaces= generate_stairs([0,0,0], 0.3, 0.1 , n_surfs); 
        case 'flat'
            surfaces= generate_contacts([0,0,0], 0.3, n_surfs);
        case 'flat_sin'
            surfaces= generate_contacts_sin([0,0,0], 0.3, n_surfs);
        case 'inclined'
            surfaces= generate_contacts([0,0,0.09], 0.3, n_surfs);
        case 'bubble'
            surfaces= generate_bubbles([0,0,0], 0.3, n_surfs);
        case 'hands'
            surfaces= generate_hands([0,0,0], 0.07, n_surfs);
        case 'box'
            surfaces= generate_box([0,0,0], 0.55, 0.4);
        case 'hands_tilt'
            surfaces= generate_hands_tilt([0,0,0], 0.15, n_surfs);

    end
    surfs=surfaces;
end
 
       

%% auxiliary functions
function surfaces = generate_contacts(start, step_lenght, n_surfaces) 
    mu          = 0.65;
    foot_lenght = 0.24;
    foot_width  = 0.18;
    feet_dist   = 0.40;

    Xr=start(1) + feet_dist/2 ;
    Xl=start(1) - feet_dist/2 ;
    Yu=start(2) + foot_lenght/2;
    Yd=start(2) - foot_lenght/2;
    z =start(3);
    r_foot = [Xr, Xr+foot_width, Xr+foot_width, Xr  ;
              Yd,            Yd,            Yu, Yu  ;
               0,             z,             z,  0 ];
    l_foot = [Xl, Xl-foot_width, Xl-foot_width, Xl  ;
              Yd,            Yd,            Yu, Yu  ;
               0,             z,             z,  0 ];
    surface1.contact_pts=r_foot;
    surface1.X=foot_width/2;
    surface1.Y=foot_lenght/2;
    surface1.mu = mu;
    surface2.contact_pts=l_foot; 
    surface2.X=foot_width/2;
    surface2.Y=foot_lenght/2;
    surface2.mu = mu;
    surfaces=[surface1,surface2];

    for step = 3 : n_surfaces 
        present = surfaces(step-1);
        new=surfaces(step-2);
        if step==n_surfaces
            new.contact_pts(2,:)=present.contact_pts(2,:);
        else
            new.contact_pts(2,:)=present.contact_pts(2,:)+ones(1,4)*step_lenght;
        end
        new.X = foot_width/2;
        new.Y = foot_lenght/2;
        new.mu = mu;
        surfaces(step)=new;
    end
end

% sinusoidal walk
function surfaces_sin = generate_contacts_sin(start, step_lenght, n_surfaces) 
    mu          = 0.65;
    foot_lenght = 0.24;
    foot_width  = 0.18;
    feet_dist   = 0.40;

    Xr=start(1) + feet_dist/2 ;
    Xl=start(1) - feet_dist/2 ;
    Yu=start(2) + foot_lenght/2;
    Yd=start(2) - foot_lenght/2;
    z =start(3);
    Rz = R_z(-atan(cos(start(2))));

    r_foot = [Xr, Xr+foot_width, Xr+foot_width, Xr  ;
              Yd,            Yd,            Yu, Yu  ;
               0,             z,             z,  0 ];
    l_foot = [Xl, Xl-foot_width, Xl-foot_width, Xl  ;
              Yd,            Yd,            Yu, Yu  ;
               0,             z,             z,  0 ];
    surface1.contact_pts=r_foot;
    surface1.X=foot_width/2;
    surface1.Y=foot_lenght/2;
    surface1.mu = mu;
    surface2.contact_pts=l_foot; 
    surface2.X=foot_width/2;
    surface2.Y=foot_lenght/2;
    surface2.mu = mu;
    surfaces=[surface1,surface2];

    surface1_sin = surface1;
    surface1_sin.contact_pts = Rz*surface1_sin.contact_pts;
    surface2_sin = surface2;
    surface2_sin.contact_pts = Rz*surface2_sin.contact_pts;
    surfaces_sin = [surface1_sin, surface2_sin];

    for step = 3 : n_surfaces 
        present = surfaces(step-1);
        new=surfaces(step-2);
        if step==n_surfaces
            new.contact_pts(2,:)=present.contact_pts(2,:);
        else
            new.contact_pts(2,:)=present.contact_pts(2,:)+ones(1,4)*step_lenght*cos(atan(cos(mean(present.contact_pts(2,:)))));
        end
        new.X = foot_width/2;
        new.Y = foot_lenght/2;
        new.mu = mu;
        surfaces(step)=new;

        new_sin = new;
        new_sin.contact_pts(1,:)=new.contact_pts(1,:)+ones(1,4)*sin(mean(new.contact_pts(2,:)))/2;
        delta = [mean(new_sin.contact_pts(1,:)); mean(new_sin.contact_pts(2,:)); mean(new_sin.contact_pts(3,:))];

        Rz = R_z(-atan(cos(mean(new.contact_pts(2,:)))));
        new_sin.contact_pts = Rz*(new_sin.contact_pts-[delta,delta,delta,delta])+[delta,delta,delta,delta];
        
        surfaces_sin(step) = new_sin;
    end
end

% stairs
function surfaces = generate_stairs(start, step_lenght, z_stair, n_surfaces) 
    mu          = 0.65;
    foot_lenght = 0.24;
    foot_width  = 0.18;
    feet_dist   = 0.40;

    Xr=start(1) + feet_dist/2 ;
    Xl=start(1) - feet_dist/2 ;
    Yu=start(2) + foot_lenght/2;
    Yd=start(2) - foot_lenght/2;
    z =start(3);
    r_foot = [Xr, Xr+foot_width, Xr+foot_width, Xr  ;
              Yd,            Yd,            Yu, Yu  ;
               z,             z,             z,  z ];
    l_foot = [Xl, Xl-foot_width, Xl-foot_width, Xl  ;
              Yd,            Yd,            Yu, Yu  ;
               z,             z,             z,  z ];
    surface1.contact_pts=r_foot;
    surface1.X=foot_width/2;
    surface1.Y=foot_lenght/2;
    surface1.mu = mu;
    surface2.contact_pts=l_foot; 
    surface2.X=foot_width/2;
    surface2.Y=foot_lenght/2;
    surface2.mu = mu;
    surfaces=[surface1,surface2];

    for step = 3 : n_surfaces 
        present = surfaces(step-1);
        new=surfaces(step-2);
        if step==n_surfaces
            new.contact_pts(2,:)=present.contact_pts(2,:);
            new.contact_pts(3,:)=present.contact_pts(3,:);
        else
            new.contact_pts(2,:)=present.contact_pts(2,:)+ones(1,4)*step_lenght;
            new.contact_pts(3,:) = present.contact_pts(3,:) + z_stair*ones(1,4);
        end
        new.X = foot_width/2;
        new.Y = foot_lenght/2;
        new.mu = mu;
        surfaces(step)=new;
    end
end

%bubbles
function surfaces = generate_bubbles(start, step_lenght, n_surfaces)
    syms x y z_sym dz_sym real
    z_sym = 0.4*(cos(2*pi*y/3.6))-0.4;

    dz_sym = jacobian(z_sym, [x,y]); 
    
    
    mu          = 0.65;
    foot_lenght = 0.24;
    foot_width  = 0.18;
    feet_dist   = 0.40;

    Xr=start(1) + feet_dist/2 ;
    Xl=start(1) - feet_dist/2 ;
    Yu=start(2) + foot_lenght/2;
    Yd=start(2) - foot_lenght/2;
    
    r_foot = [Xr, Xr+foot_width, Xr+foot_width, Xr  ;
              Yd,            Yd,            Yu, Yu  ;
               0,             0,             0,  0 ];
    l_foot = [Xl, Xl-foot_width, Xl-foot_width, Xl  ;
              Yd,            Yd,            Yu, Yu  ;
               0,             0,             0,  0 ];
    x=mean(r_foot(1,:));
    y=mean(r_foot(2,:));
    z_m = eval(z_sym);
    dz=eval(dz_sym);
    z = z_m + (r_foot(1,:)-ones(1,4)*x)*dz(1) + (r_foot(2,:)-ones(1,4)*y)*dz(2);
    r_foot(3,:)=z;

    x=mean(l_foot(1,:));
    y=mean(l_foot(2,:));
    z_m = eval(z_sym);
    dz=eval(dz_sym);
    z = z_m + (l_foot(1,:)-ones(1,4)*x)*dz(1) + (l_foot(2,:)-ones(1,4)*y)*dz(2);
    l_foot(3,:)=z;

    surface1.contact_pts=r_foot;
    surface1.X=foot_width/2;
    surface1.Y=foot_lenght/2;
    surface1.mu = mu;
    surface2.contact_pts=l_foot; 
    surface2.X=foot_width/2;
    surface2.Y=foot_lenght/2;
    surface2.mu = mu;
    surfaces=[surface1,surface2];

    for step = 3 : n_surfaces 
        present = surfaces(step-1);
        new=surfaces(step-2);
        if step==n_surfaces
            new.contact_pts(2,:)=present.contact_pts(2,:);
        else
            new.contact_pts(2,:)=present.contact_pts(2,:)+ones(1,4)*step_lenght;
        end

        x=mean(new.contact_pts(1,:));
        y=mean(new.contact_pts(2,:));
        z_m = eval(z_sym);
        dz=eval(dz_sym);
        z = z_m + (new.contact_pts(1,:)-ones(1,4)*x)*dz(1) + (new.contact_pts(2,:)-ones(1,4)*y)*dz(2);
        new.contact_pts(3,:)=z;

        new.X = foot_width/2;
        new.Y = foot_lenght/2;
        new.mu = mu;
        surfaces(step)=new;
    end


    [x_t,y_t] = meshgrid(-1:0.1:1, 0:0.1:4); 
    z_t = 0.4*(cos(2*pi*y_t/3.6))-0.4;
    mes = mesh(x_t,y_t,z_t, 'FaceAlpha',0.2);
    mes.EdgeColor = [0.3010 0.7450 0.9330];
    mes.FaceColor = [0.3010 0.7450 0.9330];
    hold on
end

% hands
function surfaces = generate_hands(start, step_lenght, n_surfaces) 
    
    mu          = 0.65;
    foot_lenght = 0.24;
    foot_width  = 0.18;
    hand_lenght = 0.18;
    hand_width  = 0.18;
    feet_dist   = 0.40;
    tilt        = 0;
    wall_dist   = 0.30;
    hand_feet_dist = 1.0;

    Xr=start(1) + wall_dist/2 ;
    Xl=start(1) - wall_dist/2 ;
    Yfu=start(2) + foot_lenght/2;
    Yfd=start(2) - foot_lenght/2;
    Yhu=start(2) + hand_lenght/2;
    Yhd=start(2) - hand_lenght/2;
    Zh=start(3) + hand_feet_dist/2 ;
    Zf=start(3) - hand_feet_dist/2 ;
    
    r_foot = [       Xr-tilt,        Xr+tilt,        Xr+tilt,       Xr-tilt  ;
                         Yfd,            Yfd,            Yfu,           Yfu  ;
               Zf-foot_width,             Zf,             Zf, Zf-foot_width ];
    l_foot = [       Xl-tilt,        Xl+tilt,        Xl+tilt,        Xl-tilt ;
                         Yfd,            Yfd,            Yfu,            Yfu ;
                          Zf,  Zf-foot_width,  Zf-foot_width,             Zf];

    r_hand = [       Xr-tilt,        Xr+tilt,        Xr+tilt,        Xr-tilt ;
                         Yhd,            Yhd,            Yhu,            Yhu ;
                          Zh,  Zh+hand_width,  Zh+hand_width,             Zh];
    l_hand = [       Xl-tilt,        Xl+tilt,        Xl+tilt,        Xl-tilt ;
                         Yhd,            Yhd,            Yhu,            Yhu ;
               Zh+hand_width,             Zh,             Zh,  Zh+hand_width];

    surface1.contact_pts=r_foot;
    surface1.X=foot_width/2;
    surface1.Y=foot_lenght/2;
    surface1.mu = mu;
    surface2.contact_pts=l_foot; 
    surface2.X=foot_width/2;
    surface2.Y=foot_lenght/2;
    surface2.mu = mu;
    surface3.contact_pts=r_hand; 
    surface3.X=hand_width/2;
    surface3.Y=hand_lenght/2;
    surface3.mu = mu;
    surface4.contact_pts=l_hand; 
    surface4.X=hand_width/2;
    surface4.Y=hand_lenght/2;
    surface4.mu = mu;
    surfaces=[surface1,surface2,surface3,surface4];

    for step = 5 : n_surfaces 
        new=surfaces(step-4);
        new.contact_pts(3,:) = new.contact_pts(3,:) + step_lenght*ones(1,4);
        surfaces(step)=new;
    end
end

% hands tilted
function surfaces = generate_hands_tilt(start, step_lenght, n_surfaces) 
    
    mu          = 0.65;
    foot_lenght = 0.24;
    foot_width  = 0.18;
    hand_lenght = 0.18;
    hand_width  = 0.18;
    feet_dist   = 0.40;
    tilt        = 0.07;
    wall_dist   = 0.70;
    hand_feet_dist = 1.0;

    Xr=start(1) + wall_dist/2 ;
    Xl=start(1) - wall_dist/2 ;
    Yfu=start(2) + foot_lenght/2;
    Yfd=start(2) - foot_lenght/2;
    Yhu=start(2) + hand_lenght/2;
    Yhd=start(2) - hand_lenght/2;
    Zh=start(3) + hand_feet_dist/2 ;
    Zf=start(3) - hand_feet_dist/2 ;
    
    r_foot = [       Xr-tilt,        Xr+tilt,        Xr+tilt,       Xr-tilt  ;
                         Yfd,            Yfd,            Yfu,           Yfu  ;
               Zf-foot_width,             Zf,             Zf, Zf-foot_width ];
    l_foot = [       Xl-tilt,        Xl+tilt,        Xl+tilt,        Xl-tilt ;
                         Yfd,            Yfd,            Yfu,            Yfu ;
                          Zf,  Zf-foot_width,  Zf-foot_width,             Zf];

    r_hand = [       Xr-tilt,        Xr+tilt,        Xr+tilt,        Xr-tilt ;
                         Yhd,            Yhd,            Yhu,            Yhu ;
                          Zh,  Zh+hand_width,  Zh+hand_width,             Zh];
    l_hand = [       Xl-tilt,        Xl+tilt,        Xl+tilt,        Xl-tilt ;
                         Yhd,            Yhd,            Yhu,            Yhu ;
               Zh+hand_width,             Zh,             Zh,  Zh+hand_width];

    surface1.contact_pts=r_foot;
    surface1.X=foot_width/2;
    surface1.Y=foot_lenght/2;
    surface1.mu = mu;
    surface2.contact_pts=l_foot; 
    surface2.X=foot_width/2;
    surface2.Y=foot_lenght/2;
    surface2.mu = mu;
    surface3.contact_pts=r_hand; 
    surface3.X=hand_width/2;
    surface3.Y=hand_lenght/2;
    surface3.mu = mu;
    surface4.contact_pts=l_hand; 
    surface4.X=hand_width/2;
    surface4.Y=hand_lenght/2;
    surface4.mu = mu;
    surfaces=[surface1,surface2,surface3,surface4];

    for step = 5 : n_surfaces 
        new=surfaces(step-4);
        new.contact_pts(3,:) = new.contact_pts(3,:) + step_lenght*ones(1,4);
        new.contact_pts(2,:) = new.contact_pts(2,:) + step_lenght*ones(1,4);
        surfaces(step)=new;
    end
end

% box
function surfaces = generate_box(start, box_h, box_dist)
    mu          = 0.65;
    %feet parmeters
    foot_lenght = 0.24;
    foot_width  = 0.18;
    feet_dist   = 0.30;

    Xfr=start(1) + feet_dist/2 ;
    Xfl=start(1) - feet_dist/2 ;
    Yfu=start(2) + foot_lenght/2;
    Yfd=start(2) - foot_lenght/2;
    zf =start(3);

    %hands parameters
    hand_width  = 0.1;  
    hand_lenght = 0.18;
    hand_dist   = 0.40;
    
    Xhr=start(1) + hand_dist/2 ;
    Xhl=start(1) - hand_dist/2 ;
    Yhu=start(2) + hand_lenght/2 + box_dist+0.2;
    Yhd=start(2) - hand_lenght/2 + box_dist+0.2;
    zh =start(3) + box_h;

    
    
    %define feet initial position
    r_foot = [Xfr, Xfr+foot_width, Xfr+foot_width, Xfr  ;
              Yfd,            Yfd,            Yfu, Yfu  ;
                zf,             zf,             zf,   zf ];
    l_foot = [Xfl, Xfl-foot_width, Xfl-foot_width, Xfl  ;
              Yfd,            Yfd,            Yfu, Yfu  ;
                zf,             zf,             zf,   zf ];
    %define hands position
    r_hand = [Xhr, Xhr+hand_width, Xhr+hand_width, Xhr ;
              Yhd,            Yhd,            Yhu, Yhu ;
                zh,             zh,             zh,  zh ];


    %define feet final position
    r_foot_f = r_foot;
    r_foot_f(2,:)= r_foot(2,:)+ ones(1,4)*box_dist;
    r_foot_f(3,:)= r_foot(3,:)+ ones(1,4)*box_h;
    l_foot_f = l_foot;
    l_foot_f(2,:)= l_foot(2,:)+ ones(1,4)*box_dist;
    l_foot_f(3,:)= l_foot(3,:)+ ones(1,4)*box_h;

    surface1.contact_pts=r_foot;
    surface1.X=foot_width/2;
    surface1.Y=foot_lenght/2;
    surface1.mu = mu;
    surface2.contact_pts=l_foot; 
    surface2.X=foot_width/2;
    surface2.Y=foot_lenght/2;
    surface2.mu = mu;
    
    surface3.contact_pts=r_hand;
    surface3.X=hand_width/2;
    surface3.Y=hand_lenght/2;
    surface3.mu = mu;

    surface5 = surface1;
    surface5.contact_pts=r_foot_f;
    
    surface6 = surface2;
    surface6.contact_pts=l_foot_f;

    surfaces=[surface2,surface1,surface3,surface6,surface5];

end




function Rz=R_z(theta)
    Rz = [ cos(theta), -sin(theta), 0;
           sin(theta),  cos(theta), 0;
                   0 ,           0, 1];
end
