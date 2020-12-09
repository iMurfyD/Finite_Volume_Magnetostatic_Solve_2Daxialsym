# Finite_Volume_Magnetostatic_Solve_2Daxialsym

Truth simulation of magnetic cohesion between paramagnetic grains. Meant to replicate results [here](https://link.aps.org/doi/10.1103/PhysRevE.89.043306). Does an axial, radial 2D finite volume simulation of a 3D azimuthally symmetric geometry. Grains are modelled as a smoothed continuous magnetic permeability profile. 

The current bottleneck in the code is the nested for loop which creates the finite volume system (`setup_system()` function). This could definitely be improved upon with vectorization or cutting out redundant conditions (`self_to_up` and `self_to_down` should be identical).

To use this, first run `test_cases/convergence_plot.m` to figure out what grid size to run at. 

Individual test cases can be run with the `solver/calc_f_two_grain()` function. An example of doing so is given in `example_force_calculation.m`. Adjust the variable `n` to a converged grid size from before. `calc_f_FV_parallel.m` can run a range of parameters (right now set to separation distance) in parallel using the MATLAB parallel toolbox. The number of cores used needs to be set up ahead of time. 
