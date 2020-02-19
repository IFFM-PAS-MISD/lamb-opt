% copy figures to paper folder and rename

load project_paths projectroot src_path;

figs_source_folder=[projectroot,'reports',filesep,'figures',filesep];
paper_folder = 'SPIE2020';
fig_destination=[projectroot,'reports',filesep,'conference_papers',filesep,paper_folder,filesep,'figs',filesep];

modelname='SASE2_plain_weave'; % rhom
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

modelname='SASE3_plain_weave'; % rhof
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

modelname='SASE4_plain_weave'; % em
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

modelname='SASE5_plain_weave'; % ef
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

modelname='SASE6_plain_weave'; % nim
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

modelname='SASE7_plain_weave'; % nif
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

modelname='SASE8_plain_weave'; % vol
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

modelname='SASE20_plain_weave'; % C11
% figure 9a
figname=[modelname,'_angle_0_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure9a.png'],'f');

% figure 9b
figname=[modelname,'_angle_30_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure9b.png'],'f');

% figure 9c
figname=[modelname,'_angle_60_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure9c.png'],'f');

% figure 9d
figname=[modelname,'_angle_90_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure9d.png'],'f');

modelname='SASE19_plain_weave'; %C12
% figure 10a
figname=[modelname,'_angle_0_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure10a.png'],'f');

% figure 10b
figname=[modelname,'_angle_30_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure10b.png'],'f');

% figure 10c
figname=[modelname,'_angle_60_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure10c.png'],'f');

% figure 10d
figname=[modelname,'_angle_90_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure10d.png'],'f');

modelname='SASE18_plain_weave'; %C13
% figure 11a
figname=[modelname,'_angle_0_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure11a.png'],'f');

% figure 11b
figname=[modelname,'_angle_30_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure11b.png'],'f');

% figure 11c
figname=[modelname,'_angle_60_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure11c.png'],'f');

% figure 11d
figname=[modelname,'_angle_90_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure11d.png'],'f');

modelname='SASE17_plain_weave'; %C22
% figure 12a
figname=[modelname,'_angle_0_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure12a.png'],'f');

% figure 12b
figname=[modelname,'_angle_30_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure12b.png'],'f');

% figure 12c
figname=[modelname,'_angle_60_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure12c.png'],'f');

% figure 12d
figname=[modelname,'_angle_90_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure12d.png'],'f');

modelname='SASE16_plain_weave'; %C23
% figure 13a
figname=[modelname,'_angle_0_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure13a.png'],'f');

% figure 13b
figname=[modelname,'_angle_30_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure13b.png'],'f');

% figure 13c
figname=[modelname,'_angle_60_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure13c.png'],'f');

% figure 13d
figname=[modelname,'_angle_90_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure13d.png'],'f');

modelname='SASE15_plain_weave'; %C33
% figure 15a
figname=[modelname,'_angle_0_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure14a.png'],'f');

% figure 14b
figname=[modelname,'_angle_30_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure14b.png'],'f');

% figure 14c
figname=[modelname,'_angle_60_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure14c.png'],'f');

% figure 14d
figname=[modelname,'_angle_90_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure14d.png'],'f');

modelname='SASE14_plain_weave'; %C44
% figure 15a
figname=[modelname,'_angle_0_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure15a.png'],'f');

% figure 15b
figname=[modelname,'_angle_30_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure15b.png'],'f');

% figure 15c
figname=[modelname,'_angle_60_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure15c.png'],'f');

% figure 15d
figname=[modelname,'_angle_90_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure15d.png'],'f');

modelname='SASE13_plain_weave'; %C55
% figure 16a
figname=[modelname,'_angle_0_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure16a.png'],'f');

% figure 16b
figname=[modelname,'_angle_30_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure16b.png'],'f');

% figure 16c
figname=[modelname,'_angle_60_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure16c.png'],'f');

% figure 16d
figname=[modelname,'_angle_90_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure16d.png'],'f');

modelname='SASE12_plain_weave'; %C66
% figure 17a
figname=[modelname,'_angle_0_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure17a.png'],'f');

% figure 17b
figname=[modelname,'_angle_30_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure17b.png'],'f');

% figure 17c
figname=[modelname,'_angle_60_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure17c.png'],'f');

% figure 17d
figname=[modelname,'_angle_90_param_dispersion_curves_color.png'];
fig_source=[figs_source_folder,'SASE',filesep,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure17d.png'],'f');
