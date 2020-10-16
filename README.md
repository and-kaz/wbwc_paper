The repository containts scripts used for the following analyses in the accompanying preprint/paper:

### 1. Calculations of dielectric permittivity using Kirkwood-Froehlich equations:
	dielectric/
		calc_dipmoment_vol.tcl 
		Tcl script for VMD to parse MD trajectories and calculate volume density maps and dipole moment for a fixed-atom-number selection;

		example_job_input.sh
		An example of bash script to setup the actual calculations using calc_dipmoment_vol.tcl
		
		calc_eps_kff.py
		Python script to calculate dielectric permittivity from volume and dipole moment trajectories;

### 2. Kinetic modeling:
	kinetic/
		ode/
			decoding_ode_server.py
			Python script to numerically integrate ordinary differential equations of the kinetic models; 
			It accepts differential equations defining the models as separate text files (provided);
			The script also contains primitive GUI to setup and analyze the systems.

		interactive_solutions/
			contains Jupyter notebook that can interactively solve and plot the analytical solutions of the kinetic models;
			
Link to the interactive_solutions on Binder:
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/and-kaz/wbwc_paper/main?filepath=%2Fkinetic%2Finteractive_solutions%2Fwbwc_decoding_analytical_solutions.ipynb)

### 3. Sripts to setup, analyze and visualize QM/MM umbrella sampling (US) simulations:
	umbrella_sampling/
		The folder contains only the scripts for (US) calculations:
		perform alignemnt to obtain frames and images for pathCV (2 Tcl scripts for VMD);
		prepare inputs and analyze results (2 Python scripts);
		an example Slurm job script to run QM/MM US calculations on a cluster.
		It is suggested to merge this folder with the folder on Figshare:
		https://doi.org/10.6084/m9.figshare.c.5178074.v1
		which contains the files necessary to setup the calculations. It also contains final colvar trajectories for one system as an example.
