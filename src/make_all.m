% make all calculations

% data processing - transformation of wavenumber field
% from cartesian to polar coordinates at selected angles
% processing of all files in interm exp folder
CartesianToPolar;

% semi analytical spectral element (SASE) model dispersion curves 
% grid search over volume fractions and Young modulus of fibres
SASE1;

% semi analytical spectral element (SASE) model dispersion curves 
% parametric search over matrix density
SASE2;

% call objective_fun for test cases calculated in SASE1,SASE2
objective_fun_script;

% plot selected results and make figures
plot_num_exp_dispersion;

