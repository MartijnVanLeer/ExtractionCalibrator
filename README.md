Source code for 'Investigating aquitard heterogeneity by inverse groundwater modelling at a drinking water production site', Van Leer et al., (2025) https://doi.org/10.1016/j.envsoft.2025.106554

The snakemake workflow creates a MODFLOW model around a well field based on a hydrogeological model, which is then calibrated.
Heterogeneous realizations are created based on borehole data and core scale hydraulic conductivity measurements. These realizations are upscaled to the model grid.
These realizations are run in the model, and the fit to drawdowns is examined. For the best fitting models MODPATH simulations are performed. 

The file 'snakefile' shows how the scripts are connected and in which order they are performed. 
In the current state it is not straightforward to directly apply the workflow to other locations, due to the large number of site specific datasets and assumptions. 
