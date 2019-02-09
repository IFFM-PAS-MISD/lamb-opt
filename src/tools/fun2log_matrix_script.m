% convert function to logical matrix

clear all; close all;

model='SASE';
% load projectroot path
load project_paths projectroot src_path;

% create path to the numerical raw data folder
raw_data_path=fullfile( projectroot, 'data','raw','num',[model,'-out'], filesep );

% create path to the experimental interim data folder
interim_data_path=fullfile( projectroot, 'data','interim','exp', filesep );

% filename of data to be processed
%filename = 'dispersion3D-3nodes-A0-4-2MHz.mat';
%filename = 'dispersion3D-3nodes-A0-4-2MHz-v3.mat';
filename = 'dispersion3D-3nodes-A0-4-2MHz-v4.mat';

% filename of interim experimental data
filename_exp = 'Krzywa.mat';
% load raw numerical data file
raw_data_num = load([raw_data_path,filename]);
load([interim_data_path,filename_exp]);
%% START DATA PROCESSING
%%
figure('Name','Wavenumber - Frequency')
%set(gcf,'Color','w','Position',[1921 100 1280 900]);
set(gcf,'Color','w');
imagesc(f_vec/1e3, k_vec_interp, abs(Wavenumber_Frequency));   set(gca,'YDir','normal'); 

set(gca,'Fontsize',14)
xlabel('f [kHz]')
ylabel('k [1/m]')
colormap(jet(255)); axis tight; 
caxis([0 max(caxis)/50]); 
%%
f_vec_max= f_vec(end);
k_vec_interp_max = k_vec_interp(end);
[km,fn] = size(Wavenumber_Frequency);
fmin=f_vec(1);
fmax=f_vec(end);
fvec2=fmin:(fmax-fmin)/(fn-1):fmax;
% numerical dispersion curve
freq=raw_data_num.om/(2*pi);
wavenumber=raw_data_num.wavenumber;
% interpolation would be unnecessary if dispersion curves will be evaluated
% at the same frequencies as in the experiment
kvec2 = interp1(freq, wavenumber,fvec2,'linear','extrap');

hold on;
%plot(freq(1:25000)/1e3,wavenumber(1:25000),'y');
plot(fvec2/1e3,kvec2,'y');

%% Function to logical matrix
% convert numerical dispersion curve into logical matrix
%H=logical(sparse(km,fn));
H=logical(zeros(km,fn));
for i=1:fn
    [~,I] = min(abs( kvec2(i) - k_vec_interp ));
    H(I,i) = 1;
end

% apply logical matrix to experimental dispersion curve matrix
dispersion = H.*abs(Wavenumber_Frequency)/max(max(abs(Wavenumber_Frequency)));
% H = flipud(H); % only for checking figure; comment later on for algorithm
% figure;
% spy(H)
figure('Name','Wavenumber - Frequency')
%set(gcf,'Color','w','Position',[1921 100 1280 900]);
set(gcf,'Color','w');
imagesc([0 f_vec_max/1e3], [0 k_vec_interp_max], dispersion);   set(gca,'YDir','normal'); 

set(gca,'Fontsize',14)
xlabel('f [kHz]')
ylabel('k [1/m]')
colormap(jet(255)); axis tight; 
caxis([0 max(caxis)/50]); 

%
[i,j] = find(H==true);
obj=0;
for k =1:length(i)
obj = obj + (sqrt(dispersion(i(k),j(k)).^2));
end
obj
return;

%% END DATA PROCESSING

% create path to the experimental interim data folder
interim_data_path=fullfile( projectroot, 'data','interim','exp', filesep );

% filename of processed data
interim_filename = 'interim_data_testfile.txt';

% save processed data to interim (intermidiate) data folder
save([interim_data_path,interim_filename],'interim_data_exp','-ascii');

% add marker (empty file '.exist' ) indicating that 
% the processing results finished successfully and exist

output_filename_marker = [interim_data_path,'.exist'];

if ~exist(output_filename_marker, 'file' )
    fid = fopen(output_filename_marker,'w');
    fclose(fid);
end
