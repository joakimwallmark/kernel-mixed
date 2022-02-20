# kernel-mixed
For simulations using IRT to generate test data, run the files
01_real_irt_parameters.R
02_start_sim_irt.R

For simulations using splines to generate test data, the equivalent files are
01_spline_cdfs.R
02_start_sim_splines.R

# Copyright
As the raw data from the SweSAT datasets cannot be shared due to copyright reasons, mock datasets
are stored in the file mock_data.RData as a replacement to enable the entire simulation process to be mimiced.

# Replicating the simulations results presented in the article
The resulting IRT parameters from the real SweSAT datasets are stored in sim_real_item_parameters.RData. 
If replication of article simulations results are of interest, 02_start_sim_irt.R can be ran directly
by loading sim_real_item_parameters.RData.

Analogously for the spline data generating process, the splines and weights obtained using the real
datasets are stored in weights.RData and fitted_splines.RData and 02_start_sim_splines.R can be ran using
the data in these files.

# Plots and tables
For simulation result plots and tables, the code in 03_result_tables.R and 04_result_plots.R can be used for both the IRT and non-IRT simulations.
