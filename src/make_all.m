%% make all calculations

%% data processing
% 3D FFT of full wavefield data
% transformation from 2D-space-time to wavenumber-frequency
% of selected file '289x289p_HANN100_x30_10Vpp_200Hz'
spatial_to_wavenumber_wavefield_selected_script;

% 3D FFT of full wavefield data (plain weave specimen) 
% transformation from 2D-space-time to wavenumber-frequency
% of selected file '499x499p_chp200_x40_18Vpp_250Hz'
spatial_to_wavenumber_wavefield_selected_script_plain_weave;

% 3D FFT of full wavefield data (unidirectional specimen) 
% transformation from 2D-space-time to wavenumber-frequency
% of selected files '499x499p_chp200_x30_6Vpp_250Hz_100mmsv_small_uni'
% '499x499p_chp200_x40_6Vpp_250Hz_uni'
spatial_to_wavenumber_wavefield_selected_script_unidirectional;

% 3D FFT of full wavefield data (Jochen's specimen with omega stringer) 
% transformation from 2D-space-time to wavenumber-frequency
% of selected files '483x483p_CHIRP_20-500kHz_125us_6Vpp_x3_stringer_intact'
% lower left quarter
spatial_to_wavenumber_wavefield_selected_script_stringer;

% 3D FFT of full wavefield data (unidirectional composite) 
% transformation from 2D-space-time to wavenumber-frequency
% from numerical simulations
spatial_to_wavenumber_wavefield_selected_script_uni_num;

% transformation of wavenumber field
% from cartesian to polar coordinates at selected angles
% processing of all files in interm exp folder
cartesian_to_polar_wavenumber_wavefield_script;

% transformation of wavenumber field
% from cartesian to polar coordinates at selected angles
% processing of all files in interm num folder
cartesian_to_polar_wavenumber_wavefield_script_num;

% transformation of wavefield
% from cartesian to polar
% of selected file '289x289p_HANN100_x30_10Vpp_200Hz'
cartesian_to_polar_wavefield_script;

% simulated pzt signals at distance about 14 cm at angles [0:15:90]
extract_signals_from_polar_wavefield

%% modelling
% semi analytical spectral element (SASE) model dispersion curves 
% grid search over volume fractions and Young modulus of fibres
SASE1;

% semi analytical spectral element (SASE) model dispersion curves 
% parametric search over matrix density
SASE2;
SASE2_plain_weave;

% semi analytical spectral element (SASE) model dispersion curves 
% parametric search over fibre density
SASE3;
SASE3_plain_weave;

% semi analytical spectral element (SASE) model dispersion curves 
% parametric search over Young modulus of matrix
SASE4;
SASE4_plain_weave;

% semi analytical spectral element (SASE) model dispersion curves 
% parametric search over Young modulus of fibres
SASE5;
SASE5_plain_weave;

% semi analytical spectral element (SASE) model dispersion curves 
% parametric search over Poisson ratio of matrix
SASE6;
SASE6_plain_weave;

% semi analytical spectral element (SASE) model dispersion curves 
% parametric search over Poisson ratio of fibres
SASE7;
SASE7_plain_weave;

% semi analytical spectral element (SASE) model dispersion curves 
% parametric search over volume fraction of reinforcing fibres
SASE8;
SASE8_plain_weave;

% semi analytical spectral element (SASE) model dispersion curves 
% grid search over volume fractions and Young modulus of fibres
% based on simulated pzt distribution (9x9 points)
SASE9;

% semi analytical spectral element (SASE) model dispersion curves 
% grid search over volume fractions and Young modulus of fibres
% based on simulated pzt distribution (17x17 points)
SASE10;

% semi analytical spectral element (SASE) model dispersion curves 
% grid search over volume fractions and Young modulus of fibres
% plain weave fabric composite
SASE11;

% (variability range 20%)
% semi analytical spectral element (SASE) model dispersion curves 
% parametric search over C66
SASE12_plain_weave;

% semi analytical spectral element (SASE) model dispersion curves 
% parametric search over C55
SASE13_plain_weave;

% semi analytical spectral element (SASE) model dispersion curves 
% parametric search over C44
SASE14_plain_weave;

% semi analytical spectral element (SASE) model dispersion curves 
% parametric search over C33
SASE15_plain_weave;

% semi analytical spectral element (SASE) model dispersion curves 
% parametric search over C23
SASE16_plain_weave;

% semi analytical spectral element (SASE) model dispersion curves 
% parametric search over C22
SASE17_plain_weave;

% semi analytical spectral element (SASE) model dispersion curves 
% parametric search over C13
SASE18_plain_weave;

% semi analytical spectral element (SASE) model dispersion curves 
% parametric search over C12
SASE19_plain_weave;

% semi analytical spectral element (SASE) model dispersion curves 
% parametric search over C11
SASE20_plain_weave;

% semi analytical spectral element (SASE) model dispersion curves 
% parametric search over C11 (variability range 30%)
SASE21_plain_weave;

% semi analytical spectral element (SASE) model dispersion curves 
% parametric search over C12 (variability range 30%)
SASE22_plain_weave;

% semi analytical spectral element (SASE) model dispersion curves 
% parametric search over C13 (variability range 30%)
SASE23_plain_weave;

% semi analytical spectral element (SASE) model dispersion curves 
% parametric search over C22 (variability range 30%)
SASE24_plain_weave;

% semi analytical spectral element (SASE) model dispersion curves 
% parametric search over C23 (variability range 30%)
SASE25_plain_weave;

% semi analytical spectral element (SASE) model dispersion curves 
% parametric search over C33 (variability range 30%)
SASE26_plain_weave;

% semi analytical spectral element (SASE) model dispersion curves 
% parametric search over C44 (variability range 30%)
SASE27_plain_weave;

% semi analytical spectral element (SASE) model dispersion curves 
% parametric search over C55 (variability range 30%)
SASE28_plain_weave;

% semi analytical spectral element (SASE) model dispersion curves 
% parametric search over C66 (variability range 30%)
SASE29_plain_weave;

%% UNIDIRECTIONAL
% semi analytical spectral element (SASE) model dispersion curves
% parametric search over C11 (variability range 30%)
SASE30_uni;
% semi analytical spectral element (SASE) model dispersion curves
% parametric search over C12 (variability range 30%)
SASE31_uni;
% semi analytical spectral element (SASE) model dispersion curves
% parametric search over C13 (variability range 30%)
SASE32_uni;
% semi analytical spectral element (SASE) model dispersion curves
% parametric search over C22 (variability range 30%)
SASE33_uni;
% semi analytical spectral element (SASE) model dispersion curves
% parametric search over C23 (variability range 30%)
SASE34_uni;
% semi analytical spectral element (SASE) model dispersion curves
% parametric search over C33 (variability range 30%)
SASE35_uni;
% semi analytical spectral element (SASE) model dispersion curves
% parametric search over C44 (variability range 30%)
SASE36_uni;
% semi analytical spectral element (SASE) model dispersion curves
% parametric search over C55 (variability range 30%)
SASE37_uni;
% semi analytical spectral element (SASE) model dispersion curves
% parametric search over C66 (variability range 30%)
SASE38_uni;

% surface dispersion (kx-ky)
% semi analytical spectral element (SASE) model dispersion curves
% parametric search over C11 (variability range 30%)
SASE30_uni_surf;
% semi analytical spectral element (SASE) model dispersion curves
% parametric search over C12 (variability range 30%)
SASE31_uni_surf;
% semi analytical spectral element (SASE) model dispersion curves
% parametric search over C13 (variability range 30%)
SASE32_uni_surf;
% semi analytical spectral element (SASE) model dispersion curves
% parametric search over C22 (variability range 30%)
SASE33_uni_surf;
% semi analytical spectral element (SASE) model dispersion curves
% parametric search over C23 (variability range 30%)
SASE34_uni_surf;
% semi analytical spectral element (SASE) model dispersion curves
% parametric search over C33 (variability range 30%)
SASE35_uni_surf;
% semi analytical spectral element (SASE) model dispersion curves
% parametric search over C44 (variability range 30%)
SASE36_uni_surf;
% semi analytical spectral element (SASE) model dispersion curves
% parametric search over C55 (variability range 30%)
SASE37_uni_surf;
% semi analytical spectral element (SASE) model dispersion curves
% parametric search over C66 (variability range 30%)
SASE38_uni_surf;

% objective functions for full wavefield laser measurements
% call objective_fun for test cases calculated in SASE1,SASE2,SASE9,SASE10, SASE11
objective_fun_nmodes_4_script;
objective_fun_nmodes_3_script;
objective_fun_nmodes_2_script;
objective_fun_nmodes_1_script;
% objective functions for simulated pzt sensors based on laser measurements
% call objective_fun_pzt for test cases calculated in SASE1
objective_fun_pzt_nmodes_4_script;
objective_fun_pzt_nmodes_3_script;
objective_fun_pzt_nmodes_2_script;
objective_fun_pzt_nmodes_1_script;
%% genetic algorithm
% Composite_Structures_GA paper
ga_plain_weave_known_mass_50;
ga_plain_weave_C_tensor_known_mass_50;
%% for second paper 
ga_unidirectonal_C_tensor_known_mass_kx_ky;
ga_unidirectional_C_tensor_known_mass;
%% plate with stringer
ga_stringer_C_tensor_known_mass;
%% visualisation
% plot and make figure for num dispersion curves on top of experimental results
plot_num_exp_dispersion_SASE1;

% plot and make figures for dispersion curves variability depending on matrix density
plot_param_dispersion_curves_SASE2;
plot_param_dispersion_curves_SASE2_color;
plot_param_dispersion_curves_SASE2_plain_weave_color;

% plot and make figures for dispersion curves variability depending on fibre density
plot_param_dispersion_curves_SASE3;
plot_param_dispersion_curves_SASE3_color;
plot_param_dispersion_curves_SASE3_plain_weave_color;

% plot and make figures for dispersion curves variability depending on Young
% modulus of matrix
plot_param_dispersion_curves_SASE4;
plot_param_dispersion_curves_SASE4_color;
plot_param_dispersion_curves_SASE4_plain_weave_color;

% plot and make figures for dispersion curves variability depending on Young
% modulus of fibres
plot_param_dispersion_curves_SASE5;
plot_param_dispersion_curves_SASE5_color;
plot_param_dispersion_curves_SASE5_plain_weave_color;

% plot and make figures for dispersion curves variability depending on Poisson
% ratio of matrix
plot_param_dispersion_curves_SASE6;
plot_param_dispersion_curves_SASE6_color;
plot_param_dispersion_curves_SASE6_plain_weave_color;

% plot and make figures for dispersion curves variability depending on Poisson
% ratio of fibres
plot_param_dispersion_curves_SASE7;
plot_param_dispersion_curves_SASE7_color;
plot_param_dispersion_curves_SASE7_plain_weave_color;

% plot and make figures for dispersion curves variability depending on volume
% fraction
plot_param_dispersion_curves_SASE8;
plot_param_dispersion_curves_SASE8_color;
plot_param_dispersion_curves_SASE8_plain_weave_color;

% (variability range 20%)
% plot and make figures for dispersion curves variability depending on C66
plot_param_dispersion_curves_SASE12_plain_weave_color;

% plot and make figures for dispersion curves variability depending on C55
plot_param_dispersion_curves_SASE13_plain_weave_color;

% plot and make figures for dispersion curves variability depending on C44
plot_param_dispersion_curves_SASE14_plain_weave_color;

% plot and make figures for dispersion curves variability depending on C33
plot_param_dispersion_curves_SASE15_plain_weave_color;

% plot and make figures for dispersion curves variability depending on C23
plot_param_dispersion_curves_SASE16_plain_weave_color;

% plot and make figures for dispersion curves variability depending on C22
plot_param_dispersion_curves_SASE17_plain_weave_color;

% plot and make figures for dispersion curves variability depending on C13
plot_param_dispersion_curves_SASE18_plain_weave_color;

% plot and make figures for dispersion curves variability depending on C12
plot_param_dispersion_curves_SASE19_plain_weave_color;

% plot and make figures for dispersion curves variability depending on C11
plot_param_dispersion_curves_SASE20_plain_weave_color;

% (variability range 30%)
% plot and make figures for dispersion curves variability depending on C11 (variability range 30%)
plot_param_dispersion_curves_SASE21_plain_weave_color;

% plot and make figures for dispersion curves variability depending on C12 (variability range 30%)
plot_param_dispersion_curves_SASE22_plain_weave_color;

% plot and make figures for dispersion curves variability depending on C13 (variability range 30%)
plot_param_dispersion_curves_SASE23_plain_weave_color;

% plot and make figures for dispersion curves variability depending on C22 (variability range 30%)
plot_param_dispersion_curves_SASE24_plain_weave_color;

% plot and make figures for dispersion curves variability depending on C23 (variability range 30%)
plot_param_dispersion_curves_SASE25_plain_weave_color;

% plot and make figures for dispersion curves variability depending on C33 (variability range 30%)
plot_param_dispersion_curves_SASE26_plain_weave_color;

% plot and make figures for dispersion curves variability depending on C44 (variability range 30%)
plot_param_dispersion_curves_SASE27_plain_weave_color;

% plot and make figures for dispersion curves variability depending on C55 (variability range 30%)
plot_param_dispersion_curves_SASE28_plain_weave_color;

% plot and make figures for dispersion curves variability depending on C66 (variability range 30%)
plot_param_dispersion_curves_SASE29_plain_weave_color;

% UNIDIRECTIONAL
% plot and make figures for dispersion curves variability depending on C11 (variability range 30%)
plot_param_dispersion_curves_SASE30_uni_color;
% plot and make figures for dispersion curves variability depending on C12 (variability range 30%)
plot_param_dispersion_curves_SASE31_uni_color;
% plot and make figures for dispersion curves variability depending on C13 (variability range 30%)
plot_param_dispersion_curves_SASE32_uni_color;
% plot and make figures for dispersion curves variability depending on C22 (variability range 30%)
plot_param_dispersion_curves_SASE33_uni_color;
% plot and make figures for dispersion curves variability depending on C23 (variability range 30%)
plot_param_dispersion_curves_SASE34_uni_color;
% plot and make figures for dispersion curves variability depending on C33 (variability range 30%)
plot_param_dispersion_curves_SASE35_uni_color;
% plot and make figures for dispersion curves variability depending on C44 (variability range 30%)
plot_param_dispersion_curves_SASE36_uni_color;
% plot and make figures for dispersion curves variability depending on C55 (variability range 30%)
plot_param_dispersion_curves_SASE37_uni_color;
% plot and make figures for dispersion curves variability depending on C66 (variability range 30%)
plot_param_dispersion_curves_SASE38_uni_color;
% surface dispersion (kx-ky)
% plot and make figures for dispersion curves variability depending on C11 (variability range 30%)
plot_param_dispersion_surf_SASE30_uni_color;
% plot and make figures for dispersion curves variability depending on C12 (variability range 30%)
plot_param_dispersion_surf_SASE31_uni_color;
% plot and make figures for dispersion curves variability depending on C13 (variability range 30%)
plot_param_dispersion_surf_SASE32_uni_color;
% plot and make figures for dispersion curves variability depending on C22 (variability range 30%)
plot_param_dispersion_surf_SASE33_uni_color;
% plot and make figures for dispersion curves variability depending on C23 (variability range 30%)
plot_param_dispersion_surf_SASE34_uni_color;
% plot and make figures for dispersion curves variability depending on C33 (variability range 30%)
plot_param_dispersion_surf_SASE35_uni_color;
% plot and make figures for dispersion curves variability depending on C44 (variability range 30%)
plot_param_dispersion_surf_SASE36_uni_color;
% plot and make figures for dispersion curves variability depending on C55 (variability range 30%)
plot_param_dispersion_surf_SASE37_uni_color;
% plot and make figures for dispersion curves variability depending on C66 (variability range 30%)
plot_param_dispersion_surf_SASE38_uni_color;
% 
%
% plot and make figure for num dispersion curves on top of experimental results
% for simulated pzt distrubution (9x9 points)
plot_num_exp_dispersion_SASE9;

% plot and make figure for num dispersion curves on top of experimental results
% for simulated pzt distrubution (17x17 points)
plot_num_exp_dispersion_SASE10;

% plot and make figure for num dispersion curves on top of experimental results (fabric composite)
plot_num_exp_dispersion_SASE11;
% plot and make figure for optimized num dispersion curves on top of experimental results (plain weave)
plot_dispersion_curves_ga_plain_weave_known_mass;
plot_dispersion_curves_ga_plain_weave_known_mass_large;
plot_dispersion_curves_ga_plain_weave_known_mass_50;
plot_dispersion_curves_ga_plain_weave_known_mass_50_large;
plot_dispersion_curves_ga_plain_weave_C_tensor_known_mass_50;
plot_dispersion_curves_ga_plain_weave_C_tensor_known_mass_50_l;
% SH0 mode problem
%plot_dispersion_curves_ga_plain_weave_C_tensor_known_mass_SH0
plot_dispersion_curves_ga_plain_weave_C_tensor_known_mass_SH0_l;

% unidirectional
plot_dispersion_curves_ga_uni_C_tensor_known_mass_large;
plot_dispersion_curves_ga_uni_C_tensor_known_mass_small
% surface plot
plot_dispersion_surf_ga_uni_C_tensor_known_mass_large;
plot_dispersion_surf_ga_uni_C_tensor_known_mass_small;

% plot and make figure for initial num dispersion curves on top of experimental results (plain weave)
plot_dispersion_curves_ga_plain_weave_known_mass_initial;
plot_dispersion_curves_ga_plain_weave_known_mass_initial_large;

% plot experimental dispersion curves for various excitation signals
plot_exp_dispersion;
plot_exp_dispersion_plain_weave;
% plot chirp signals for reports and presentations
plot_chirp_signals;

% plot objective function scores for SASE1 comparing laser and simulated pzt
plot_objective_function_score; % vibrant colour palette
plot_objective_function_score_high_contrast; % high-contrast palette

% plot exemplary dispersion curves showing modes tracking problem
plot_exemplary_dispersion_curves_SASE2;
plot_exemplary_dispersion_curves_SASE2_sorted;

% Composite_Structures_GA paper
% tables
plot_constants_histogram_50;

% GA convergence fig (commented at the end)
%ga_plain_weave_C_tensor_known_mass.m;

% copy figs to paper folder
copy_fig_Composite_Structures_GA;
copy_fig_SPIE2020;