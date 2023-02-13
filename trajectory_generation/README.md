# Trajectory Generation scripts
Download Casadi for MATLAB: https://web.casadi.org/get/

flat_trajectory_giw.m :                                                    Generates wrench cones using the function gravito_inertial_wrench.m,
                                      solve an optimization problem using Casadi library using system dynamics and wrench cones as constraints,
              in this case we consider stances on the same flat ground and solve the problem as a single optimization for the whole trajectory. 


nonflat_trajectory_giw.m :                                                         Generates wrench cones using the function gravito_inertial_wrench.m,
                                              solve an optimization problem using Casadi library using system dynamics and wrench cones as constraints,
                 in this case we consider stances on different inclined planes and solve the problem as a single optimization for the whole trajectory. 
                 
