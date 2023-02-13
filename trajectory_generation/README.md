# Trajectory Generation scripts
Download Casadi for MATLAB: https://web.casadi.org/get/

### flat_trajectory_giw.m :                                                    
Generates wrench cones using the function gravito_inertial_wrench.m,

solve an optimization problem using Casadi library using system dynamics and wrench cones as constraints,

in this case we consider stances on the same flat ground and solve the problem as a single optimization for the whole trajectory. 


### nonflat_trajectory_giw.m :                                                         
Generates wrench cones using the function gravito_inertial_wrench.m,

solve an optimization problem using Casadi library using system dynamics and wrench cones as constraints,

in this case we consider stances on different inclined planes and solve the problem as a single optimization for the whole trajectory.

PROBLEM : We can't find a solution for more than 3 different stances

### iterative_two_stances.m :
The problem of multiple stances optimization is solved by optimizing the trajectory considering two stances at a time,

The two stances problem is solved n_stances-1 times varing initial, final positions and Gravito-inertial wrench cones.
                                        
### iterative_one_stance.m :
The problem of multiple stances optimization is solved by optimizing the trajectory considering one stance at a time,

The one stances problem is solved n_stances times varing initial, final positions and Gravito-inertial wrench cones.
       
