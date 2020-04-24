% copy figures to paper folder and rename

load project_paths projectroot src_path;

figs_source_folder=[projectroot,'reports',filesep,'figures',filesep];
paper_folder = 'Composite_Structures_GA_R2';
fig_destination=[projectroot,'reports',filesep,'journal_papers',filesep,paper_folder,filesep,'figs',filesep];

%% parametric studies

modelname='SASE21_plain_weave'; % C11
% figure 2a
figname=[modelname,'_angle_0_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure2a.png'],'f');

% figure 2b
figname=[modelname,'_angle_45_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure2b.png'],'f');

modelname='SASE22_plain_weave'; %C12
% figure 3a
figname=[modelname,'_angle_0_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure3a.png'],'f');

% figure 3b
figname=[modelname,'_angle_45_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure3b.png'],'f');

modelname='SASE23_plain_weave'; %C13
% figure 4a
figname=[modelname,'_angle_0_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure4a.png'],'f');

% figure 4b
figname=[modelname,'_angle_45_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure4b.png'],'f');

modelname='SASE24_plain_weave'; %C22
% figure 5a
figname=[modelname,'_angle_0_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure5a.png'],'f');

% figure 5b
figname=[modelname,'_angle_45_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure5b.png'],'f');


modelname='SASE25_plain_weave'; %C23
% figure 6a
figname=[modelname,'_angle_0_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure6a.png'],'f');

% figure 6b
figname=[modelname,'_angle_45_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure6b.png'],'f');

modelname='SASE26_plain_weave'; %C33
% figure 7a
figname=[modelname,'_angle_0_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure7a.png'],'f');

% figure 7b
figname=[modelname,'_angle_45_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure7b.png'],'f');


modelname='SASE27_plain_weave'; %C44
% figure 8a
figname=[modelname,'_angle_0_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure8a.png'],'f');

% figure 8b
figname=[modelname,'_angle_45_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure8b.png'],'f');

modelname='SASE28_plain_weave'; %C55
% figure 14a
figname=[modelname,'_angle_0_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure14a.png'],'f');

% figure 14b
figname=[modelname,'_angle_45_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure14b.png'],'f');

modelname='SASE29_plain_weave'; %C66
% figure 15a
figname=[modelname,'_angle_0_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure15a.png'],'f');

% figure 15b
figname=[modelname,'_angle_45_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure15b.png'],'f');

%% genetic algorithm
modelname='ga_plain_weave_known_mass';
% figure 10a
figname=[modelname,'_angle_60_dispersion_curves_initial.png'];
fig_source=[figs_source_folder,'genetic_algorithm',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure10a.png'],'f');

modelname='ga_plain_weave_known_mass_50';
% figure 10b
figname=[modelname,'_angle_60_dispersion_curves_test_case_4.png'];
fig_source=[figs_source_folder,'genetic_algorithm',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure10b.png'],'f');

% figure 11a
figname=[modelname,'_angle_0_dispersion_curves_test_case_4.png'];
fig_source=[figs_source_folder,'genetic_algorithm',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure11a.png'],'f');

% figure 11b
figname=[modelname,'_angle_15_dispersion_curves_test_case_4.png'];
fig_source=[figs_source_folder,'genetic_algorithm',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure11b.png'],'f');

% figure 11c
figname=[modelname,'_angle_30_dispersion_curves_test_case_4.png'];
fig_source=[figs_source_folder,'genetic_algorithm',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure11c.png'],'f');

% figure 11d
figname=[modelname,'_angle_45_dispersion_curves_test_case_4.png'];
fig_source=[figs_source_folder,'genetic_algorithm',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure11d.png'],'f');

% figure 11e
figname=[modelname,'_angle_75_dispersion_curves_test_case_4.png'];
fig_source=[figs_source_folder,'genetic_algorithm',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure11e.png'],'f');

% figure 11f
figname=[modelname,'_angle_90_dispersion_curves_test_case_4.png'];
fig_source=[figs_source_folder,'genetic_algorithm',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure11f.png'],'f');

modelname='ga_plain_weave_C_tensor_known_mass_50';

% figure 12a
figname=[modelname,'_angle_0_dispersion_curves_test_case_2_large.png'];
fig_source=[figs_source_folder,'genetic_algorithm',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure12a.png'],'f');

% figure 12b
figname=[modelname,'_angle_15_dispersion_curves_test_case_2_large.png'];
fig_source=[figs_source_folder,'genetic_algorithm',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure12b.png'],'f');

% figure 13a
figname=[modelname,'_angle_0_dispersion_curves_test_case_6.png'];
fig_source=[figs_source_folder,'genetic_algorithm',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure13a.png'],'f');

% figure 13b
figname=[modelname,'_angle_15_dispersion_curves_test_case_6.png'];
fig_source=[figs_source_folder,'genetic_algorithm',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure13b.png'],'f');

% figure 13c
figname=[modelname,'_angle_30_dispersion_curves_test_case_6.png'];
fig_source=[figs_source_folder,'genetic_algorithm',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure13c.png'],'f');

% figure 13d
figname=[modelname,'_angle_45_dispersion_curves_test_case_6.png'];
fig_source=[figs_source_folder,'genetic_algorithm',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure13d.png'],'f');

% figure 13e
figname=[modelname,'_angle_75_dispersion_curves_test_case_6.png'];
fig_source=[figs_source_folder,'genetic_algorithm',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure13e.png'],'f');

% figure 13f
figname=[modelname,'_angle_90_dispersion_curves_test_case_6.png'];
fig_source=[figs_source_folder,'genetic_algorithm',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure13f.png'],'f');