# Wrench Cone Generation scripts

contact_cone.m : contains a simple attempt of a single contact point on a surface.

multiple_points.m :                                extends the simple contact to a convex surface defined by 5 vertex, 
                                                                 computes the span and face form of the contact cones,
                                                               defines the stance given by the contact on the surface,
                    computes the correspondent Gravito-Inertial cone applied to the Center of pressure of the surface.

multiple_surfaces.m:           extends multiple contact considering the same pattern of contacts on multiple surfaces,
                                                                         computes span and face form for each surface,
                                                                   defines the stance considering all contact sufaces,
                                          computes Gravito_Inertial wrench cone applied to an arbitrarily choosen CoM.
gravito_inertial_wrench.m:                                    it is basically multiple_surfaces but in function shape,
                                                                                                input array of struct, 
                                                                                 output wrench cone and its face form,
                                the struct represents a surface by its contact points and static friction coefficient.     
