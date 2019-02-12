% make all calculations

% data processing

CartesianToPolar;

% semi analytical spectral element (SASE) model dispersion curves 
% grid search over volume fractions and Youn modulus of fibres

SASE1;

% call objective_fun for test cases calculated in SASE1
objective_fun_test_script;

% plot selected results and make figures
plot_num_exp_dispersion;

