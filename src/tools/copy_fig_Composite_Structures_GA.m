% copy figures to paper folder and rename

load project_paths projectroot src_path;

figs_source_folder=[projectroot,'reports',filesep,'figures',filesep];
paper_folder = 'Composite_Structures_GA';
fig_destination=[projectroot,'reports',filesep,'journal_papers',filesep,paper_folder,filesep,'figs',filesep];

modelname='SASE2_plain_weave';
% figure 2a
figname=[modelname,'_angle_0_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);

% figure 2b
figname=[modelname,'_angle_30_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);

modelname='SASE3_plain_weave';
% figure 3a
figname=[modelname,'_angle_0_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);

% figure 3b
figname=[modelname,'_angle_30_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);

modelname='SASE4_plain_weave';
% figure 4a
figname=[modelname,'_angle_0_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);

% figure 4b
figname=[modelname,'_angle_30_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);

modelname='SASE5_plain_weave';
% figure 5a
figname=[modelname,'_angle_0_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);

% figure 5b
figname=[modelname,'_angle_30_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);

modelname='SASE6_plain_weave';
% figure 6a
figname=[modelname,'_angle_0_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);

% figure 6b
figname=[modelname,'_angle_30_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);

modelname='SASE7_plain_weave';
% figure 7a
figname=[modelname,'_angle_0_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);

% figure 7b
figname=[modelname,'_angle_30_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);

modelname='SASE8_plain_weave';
% figure 8a
figname=[modelname,'_angle_0_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);

% figure 8b
figname=[modelname,'_angle_30_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);