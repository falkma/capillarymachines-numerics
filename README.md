## Microscopic manipulation with capillary machines (numerics)

The code uploaded here is for replicating the simulation results of the paper "3D-printed machines that manipulate microscopic objects using capillary forces" by Zeng*, Faaborg*, et al. Almost all code is written for MATLAB (v2018b) with the PDE Toolbox installed. A few scripts are written in Python (v3.7).

The code is divided in three sections, which we describe in more detail:

### Ratchet

- 01_ratchet_sweep_20.m: Initial script to run, generates an energy landscape across the two geometrical parameters &alpha; and &beta; of the ratchet device described in the paper. A steep contact angle of 20&deg; is assumed between the channel wall and the meniscus.
- 02_ratchet_sweep_70.m: Second script to run, generates an energy landscape across two geometrical parameters &alpha; and &beta; which describe the ratchet geometry. A shallow contact angle of 70&deg; is assumed between the channel wall and the meniscus.
- 03_navigating_landscape: Third script to run, implements a simple greedy descent procedure to track the evolution of the minimum energy &beta; as &alpha; is changed.

### Taichi

- 01_taichi_forward_sim.m: Initial script to run, simulates trajectory of floats on the upward, braiding stroke of the taichi device described in the paper.
- 02_taichi_forward_sim_refine.m: Second script to run, refines the trajectory of the upwards stroke along regions with high variability in the trajectories.
- 03_taichi_reverse_sim.m: Third script to run, simulates trajectory of floats on the downwards, reset stroke of the taichi device.
- 04_taichi_forward_viz.m: Fourth script to run, loads and visualize the trajectories of floats in the refined forward taichi simulation.
- 05_taichi_reverse_vid.m: Fifth script to run, loads and visualize the trajectories of floats in the reverse taichi simulation.
- 06_make_movie_full.m: Sixth script to run, produces movie from saved images of taichi simulation.

### Design Space

Contains three folders, one for each parameter sweep over the three design parameters we consider for the capillary trap geometry: (1) normal force on float, (2) float radius, and (3) contact slope at the channel boundary. See paper for further description of these variables. The first script in each folder runs the parameter sweep and outputs energies as a function of the considered design parameter. The second script plots these numerical results against theoretical predictions derived in the supplement to the paper.
