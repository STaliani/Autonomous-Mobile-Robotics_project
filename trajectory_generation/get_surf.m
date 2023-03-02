% Select the type of surface: four_points_flat, four_points_nonflat
function surfs = get_surf(type_of_surface, n_surfs)

    switch type_of_surface
        case 'four_points_flat'
            surfaces= generate_contacts([0,0,0], 0.6, n_surfs);
        case 'four_points_nonflat'
            surfaces= generate_contacts([0,0,0.09], 0.6, n_surfs);
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