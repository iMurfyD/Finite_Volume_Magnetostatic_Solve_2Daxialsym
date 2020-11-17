# Finite_Volume_Magnetostatic_Solve_2Daxialsym

Truth simulation of magnetic cohesion between paramagnetic grains. Meant to replicate results [here](https://link.aps.org/doi/10.1103/PhysRevE.89.043306). Does an axial, radial 2D finite volume simulation of a 3D azimuthally symmetric geometry. Grains are modelled as a smoothed continuous magnetic permeability profile. 

The numerical force integration method isn't great (noticeably more noisey than the actual Maxwell stress tensor calculation). Bear that in mind when running. 

The current bottleneck in the code is the nested for loop which creates the finite volume system (`setup_system_sparse()` function). This could definitely be improved upon with vectorization or cutting out redundant conditions (`self_to_up` and `self_to_down` should be identical).

To use this, first run `convergence_plot.m` to figure out what grid size to run at. Then `calc_truth_f_circum_method.m` can run individual parameters. `calc_f_FV_parallel.m` can run a range of parameters (right now set to separation distance) in parallel using the MATLAB parallel toolbox. The number of cores used needs to be set up ahead of time.
