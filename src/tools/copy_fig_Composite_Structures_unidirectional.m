% copy figures to paper folder and rename

load project_paths projectroot src_path;

figs_source_folder=[projectroot,'reports',filesep,'figures',filesep];
paper_folder = 'Composite_Structures_unidirectional';
fig_destination=[projectroot,'reports',filesep,'journal_papers',filesep,paper_folder,filesep,'figs',filesep];

%% 3D slices
modelfolder='Composite_Structures_uni';
foldername='plot_exp_3Dslices_Composite_Structures'; 
% figure 2
figname=['angle_slice2.png'];
fig_source=[figs_source_folder,modelfolder,filesep,foldername,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure2.png'],'f');
% figure 3
figname=['freq_slice2.png'];
fig_source=[figs_source_folder,modelfolder,filesep,foldername,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure3.png'],'f');

%% homogenised vs optimised
% homogenised
modelfolder='genetic_algorithm';
foldername='ga_uni_SHMII_2021_homogenized'; 
% figure 4a
figname=[foldername,'_angle_0_dispersion_curves_test_case_homogenized_small.png'];
fig_source=[figs_source_folder,modelfolder,filesep,foldername,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure4a.png'],'f');
% figure 4c
figname=[foldername,'_angle_45_dispersion_curves_test_case_homogenized_small.png'];
fig_source=[figs_source_folder,modelfolder,filesep,foldername,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure4c.png'],'f');
% figure 4e
figname=[foldername,'_angle_60_dispersion_curves_test_case_homogenized_small.png'];
fig_source=[figs_source_folder,modelfolder,filesep,foldername,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure4e.png'],'f');
% figure 4g
figname=[foldername,'_angle_90_dispersion_curves_test_case_homogenized_small.png'];
fig_source=[figs_source_folder,modelfolder,filesep,foldername,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure4g.png'],'f');

% optimised
modelfolder='genetic_algorithm';
foldername='ga_uni_SHMII_2021_optimized'; 
% figure 4b
figname=['ga_unidirectional_C_tensor_known_mass_mut_rnd_offspring_2lay6_angle_0_dispersion_curves_test_case_1_small.png'];
fig_source=[figs_source_folder,modelfolder,filesep,foldername,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure4b.png'],'f');
% figure 4d
figname=['ga_unidirectional_C_tensor_known_mass_mut_rnd_offspring_2lay6_angle_45_dispersion_curves_test_case_1_small.png'];
fig_source=[figs_source_folder,modelfolder,filesep,foldername,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure4d.png'],'f');
% figure 4f
figname=['ga_unidirectional_C_tensor_known_mass_mut_rnd_offspring_2lay6_angle_60_dispersion_curves_test_case_1_small.png'];
fig_source=[figs_source_folder,modelfolder,filesep,foldername,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure4f.png'],'f');
% figure 4h
figname=['ga_unidirectional_C_tensor_known_mass_mut_rnd_offspring_2lay6_angle_90_dispersion_curves_test_case_1_small.png'];
fig_source=[figs_source_folder,modelfolder,filesep,foldername,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure4h.png'],'f');
%% f-slice
% optimised
modelfolder='genetic_algorithm';
foldername='ga_unidirectional_C_tensor_known_mass_kx_ky'; 
% figure 5a
figname=[foldername,'_frequency_3_dispersion_surf_test_case_2_small.png'];
fig_source=[figs_source_folder,modelfolder,filesep,foldername,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure5a.png'],'f');
% figure 5b
figname=[foldername,'_frequency_4_dispersion_surf_test_case_2_small.png'];
fig_source=[figs_source_folder,modelfolder,filesep,foldername,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure5b.png'],'f');
% figure 5c
figname=[foldername,'_frequency_5_dispersion_surf_test_case_2_small.png'];
fig_source=[figs_source_folder,modelfolder,filesep,foldername,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure5c.png'],'f');
% figure 5d
figname=[foldername,'_frequency_7_dispersion_surf_test_case_2_small.png'];
fig_source=[figs_source_folder,modelfolder,filesep,foldername,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure5d.png'],'f');


%% SASE vs SEM
% SASE optimised SEM homogenised 
modelfolder='Composite_Structures_uni';
foldername='ga_uni_CS_2021_homogenized_num'; 
% figure 6a
figname=[foldername,'_angle_0_dispersion_curves_test_case_1_small.png'];
fig_source=[figs_source_folder,modelfolder,filesep,foldername,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure6a.png'],'f');
% figure 6c
figname=[foldername,'_angle_45_dispersion_curves_test_case_1_small.png'];
fig_source=[figs_source_folder,modelfolder,filesep,foldername,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure6c.png'],'f');
% figure 6e
figname=[foldername,'_angle_60_dispersion_curves_test_case_1_small.png'];
fig_source=[figs_source_folder,modelfolder,filesep,foldername,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure6e.png'],'f');

% SASE optimised SEM optimised
modelfolder='Composite_Structures_uni';
foldername='ga_uni_CS_2021_optimized_num'; 
% figure 6b
figname=[foldername,'_angle_0_dispersion_curves_test_case_1_small.png'];
fig_source=[figs_source_folder,modelfolder,filesep,foldername,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure6b.png'],'f');
% figure 6d
figname=[foldername,'_angle_45_dispersion_curves_test_case_1_small.png'];
fig_source=[figs_source_folder,modelfolder,filesep,foldername,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure6d.png'],'f');
% figure 6f
figname=[foldername,'_angle_60_dispersion_curves_test_case_1_small.png'];
fig_source=[figs_source_folder,modelfolder,filesep,foldername,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure6f.png'],'f');

% figure 7a
% ma-shm project
figname='Vz_20_frame129_bottom_500.00.png';
fig_source=['/home/pkudela/work/projects/nawa-bekker/ma-shm/reports/figures/flat_shell/flat_shell_SHMII2021_out/20_output/',figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure7a.png'],'f');

% figure 7b
% ma-shm project
figname='Vz_17_frame129_bottom_500.00.png';
fig_source=['/home/pkudela/work/projects/nawa-bekker/ma-shm/reports/figures/flat_shell/flat_shell_SHMII2021_out/17_output/',figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure7b.png'],'f');

% figure 7c
% ma-shm project
figname='492x492p_50kHz_5HC_x20_15Vpp_frame257_500.00.png';
fig_source=['/home/pkudela/work/projects/nawa-bekker/ma-shm/reports/figures/EWSHM2020_exp_wavefield_out/2_output/',figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure7c.png'],'f');

%% f-slice comparison with Moll OGW
% optimised
modelfolder='genetic_algorithm';
foldername='ga_stringer_C_tensor_known_mass_kx_ky_2'; 
% figure 8a
figname=[foldername,'_frequency_3_dispersion_surf_test_case_1_small2.png'];
fig_source=[figs_source_folder,modelfolder,filesep,foldername,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure8a.png'],'f');
% figure 8b
figname=[foldername,'_frequency_4_dispersion_surf_test_case_1_small2.png'];
fig_source=[figs_source_folder,modelfolder,filesep,foldername,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure8b.png'],'f');
% figure 8c
figname=[foldername,'_frequency_5_dispersion_surf_test_case_1_small2.png'];
fig_source=[figs_source_folder,modelfolder,filesep,foldername,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure8c.png'],'f');
% figure 8d
figname=[foldername,'_frequency_7_dispersion_surf_test_case_1_small2.png'];
fig_source=[figs_source_folder,modelfolder,filesep,foldername,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure8d.png'],'f');

