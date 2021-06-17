
disp('loading wavefield data');


wavefield_name='333x333p_100kHz_5HC_10Vpp_x10_pzt';
load(['/home/pkudela/work/projects/opus15/lamb-opt/data/processed/exp/A0_S0_SH0_mode_filtering_3D/',wavefield_name,'_A0']);% Data_A0, WL, time
load(['/home/pkudela/work/projects/opus15/lamb-opt/data/processed/exp/A0_S0_SH0_mode_filtering_3D/',wavefield_name,'_S0']);% Data_S0, WL, time
load(['/home/pkudela/work/projects/opus15/lamb-opt/data/processed/exp/A0_S0_SH0_mode_filtering_3D/',wavefield_name,'_SH0']);% Data_SH0, WL, time

Width=WL(1);
Length=WL(2);
%RMS_A0= squeeze(sum(Data_A0.^2,3));
beta = 0:5:90;
% central point
I=167; J=168;

[nY,nX,nT] = size(Data_A0);
dx=Length/(nX-1);
dy=Width/(nY-1);
% dimensions of the quarter
Lx=(nX-I)*dx;
Ly=(nY-J)*dy;   
[Data_polar_A0,number_of_points,radius] = cartesian_to_polar_wavefield2(Data_A0(I:end,J:end,:),Lx,Ly,beta); 
[Data_polar_S0,number_of_points,radius] = cartesian_to_polar_wavefield2(Data_S0(I:end,J:end,:),Lx,Ly,beta); 
[Data_polar_SH0,number_of_points,radius] = cartesian_to_polar_wavefield2(Data_SH0(I:end,J:end,:),Lx,Ly,beta); 
Data_polar_A0_magnitude=squeeze(max(abs(Data_polar_A0),[],3));
Data_polar_S0_magnitude=squeeze(max(abs(Data_polar_S0),[],3));
Data_polar_SH0_magnitude=squeeze(max(abs(Data_polar_SH0),[],3));
%% dispersion curves
% modelfolder = 'genetic_algorithm';
% modelname = 'ga_unidirectional_C_tensor_known_mass_kx_ky';
% % path to dispersion curves data
% dispersion_output_path = prepare_figure_paths(modelfolder,modelname);
% disp('loading dispersion curves');
% load([dispersion_output_path,'dispersion3D_uni1']);% 'wavenumber1','CG1','FREQ1' - dispersion curves for optimized constants
nsT = 1024;
Fs =  1/(time(3)-time(2));                % sampling frequency
f_vec = Fs/2*linspace(0,1,nsT/2);         % frequency vector
modelname2 = 'A0_S0_SH0_mode_filtering_3D_chirp';
output_path = prepare_folder_paths('processed','exp',modelname2);
no_of_angles=512;
beta2 = 0:90/(no_of_angles-1):90;
load([output_path,'wavenumber1q_A0']);%[wavenumber, angle]
load([output_path,'wavenumber1q_S0']);
load([output_path,'wavenumber1q_SH0']);
wavenumber1q_A0_new=interp1(beta2',permute(wavenumber1q_A0,[2 1]),beta','spline'); % interpolate on lower number of angles
wavenumber1q_S0_new=interp1(beta2',permute(wavenumber1q_S0,[2 1]),beta','spline'); % interpolate on lower number of angles
wavenumber1q_SH0_new=interp1(beta2',permute(wavenumber1q_SH0,[2 1]),beta','spline'); % interpolate on lower number of angles
wavenumber1q_A0_new = wavenumber1q_A0_new/(2*pi); % change units to [1/m]
wavenumber1q_S0_new = wavenumber1q_S0_new/(2*pi); % change units to [1/m]
wavenumber1q_SH0_new = wavenumber1q_SH0_new/(2*pi); % change units to [1/m]
% index for 100 kHz
[a,I]=min(abs(f_vec-100000));
cp_A0=f_vec(I)./wavenumber1q_A0_new(:,I); 
cp_S0=f_vec(I)./wavenumber1q_S0_new(:,I); 
cp_SH0=f_vec(I)./wavenumber1q_SH0_new(:,I); 

%% fit curves
R=linspace(0,radius,number_of_points);
% fit curves in a distance range 5 cm - 20 cm
% avoiding near field and edge reflections
% R(35) - R(134)
% grid search optimization
Rstart=35;
Rend=134;
alpha=19132;
% A0 mode
A0=[0.0005:0.00005:0.004]'; % amplitude range for optimization
f_A0=zeros(length(beta),number_of_points);
A_A0=zeros(length(beta),1);
for k=1:length(beta)
    eta0=alpha/(2*cp_A0(k));
    f0=A0./sqrt(R).*exp(-eta0.*R);
    f2=sum((f0(:,Rstart:Rend)-Data_polar_A0_magnitude(k,Rstart:Rend)).^2,2);
    [~,J] = min(f2);
    f_A0(k,:)=f0(J,:);
    A_A0(k)=A0(J);
end

% S0 mode
alpha2=1*alpha;
S0=[0.00005:0.00001:0.0009]'; % amplitude range for optimization
f_S0=zeros(length(beta),number_of_points);
A_S0=zeros(length(beta),1);
for k=1:length(beta)
    eta0=alpha/(2*cp_S0(k));
    f0=S0./sqrt(R).*exp(-eta0.*R);
    f2=sum((f0(:,Rstart:Rend)-Data_polar_S0_magnitude(k,Rstart:Rend)).^2,2);
    [~,J] = min(f2);
    f_S0(k,:)=f0(J,:);
    A_S0(k)=S0(J);
end

% SH0 mode
SH0=[0.00005:0.00001:0.001]'; % amplitude range for optimization
f_SH0=zeros(length(beta),number_of_points);
A_SH0=zeros(length(beta),1);
for k=1:length(beta)
    eta0=alpha/(2*cp_SH0(k));
    f0=SH0./sqrt(R).*exp(-eta0.*R);
    f2=sum((f0(:,Rstart:Rend)-Data_polar_SH0_magnitude(k,Rstart:Rend)).^2,2);
    [~,J] = min(f2);
    f_SH0(k,:)=f0(J,:);
    A_SH0(k)=SH0(J);
end
%% plotting
% phase velocities polar plots
figure;
polarplot(beta*pi/180,cp_A0,'ro');hold on;
polarplot(beta*pi/180,cp_S0,'bd');
polarplot(beta*pi/180,cp_SH0,'gv');
legend('A0','S0','SH0');
thetalim([0 90]);
set(gcf,'Color','w');
set(gca,'FontName','Times');
title({'c_p [m/s]'});
% eta (attenuation) polar plots
figure;
    polarplot(beta*pi/180,alpha./(2*cp_A0),'ro');hold on;
    polarplot(beta*pi/180,alpha./(2*cp_S0),'bd');
    polarplot(beta*pi/180,alpha./(2*cp_SH0),'gv');
    legend('A0','S0','SH0');
    thetalim([0 90]);
    set(gcf,'Color','w');
    set(gca,'FontName','Times');
    title({'\eta [Np/m]'});
% A0 mode
figure;
    plot(R,Data_polar_A0_magnitude(1,:),'g'); % 0 deg
    hold on;
    plot(R,Data_polar_A0_magnitude(7,:),'k'); % 30 deg
    plot(R,Data_polar_A0_magnitude(10,:),'r'); % 45 deg
    plot(R,Data_polar_A0_magnitude(14,:),'m'); % 65 deg
    plot(R,Data_polar_A0_magnitude(19,:),'b'); % 90 deg
    xlabel('d [m]');
    ylabel('A [m/s]');
    
    title('A0 mode');

    plot(R,f_A0(1,:),'g--');
    plot(R,f_A0(7,:),'k--');
    plot(R,f_A0(10,:),'r--');
    plot(R,f_A0(14,:),'m--');
    plot(R,f_A0(19,:),'b--');
    legend([num2str(beta(1)), ' deg'], [num2str(beta(7)), ' deg'], [num2str(beta(10)), ' deg'], ...
           [num2str(beta(14)), ' deg'], [num2str(beta(19)), ' deg'],...
           ['A= ',num2str(A_A0(1)),' \eta = ',num2str(alpha/(2*cp_A0(1)))],...
           ['A= ',num2str(A_A0(7)),' \eta = ',num2str(alpha/(2*cp_A0(7)))],...
           ['A= ',num2str(A_A0(10)),' \eta = ',num2str(alpha/(2*cp_A0(10)))],...
           ['A= ',num2str(A_A0(14)),' \eta = ',num2str(alpha/(2*cp_A0(14)))],...
           ['A= ',num2str(A_A0(19)),' \eta = ',num2str(alpha/(2*cp_A0(19)))]);
    xlim([0 0.2]);ylim([0 0.03]);


% S0 mode
figure;
    plot(R,Data_polar_S0_magnitude(1,:),'g');
    hold on;
    plot(R,Data_polar_S0_magnitude(7,:),'k');
    plot(R,Data_polar_S0_magnitude(10,:),'r');
    plot(R,Data_polar_S0_magnitude(14,:),'m');
    plot(R,Data_polar_S0_magnitude(19,:),'b');


    plot(R,f_S0(1,:),'g--');
    plot(R,f_S0(7,:),'k--');
    plot(R,f_S0(10,:),'r--');
    plot(R,f_S0(14,:),'m--');
    plot(R,f_S0(19,:),'b--');

    xlabel('d [m]');
    ylabel('A [m/s]');
    legend([num2str(beta(1)), ' deg'], [num2str(beta(7)), ' deg'], [num2str(beta(10)), ' deg'], ...
           [num2str(beta(14)), ' deg'], [num2str(beta(19)), ' deg'],...
           ['A= ',num2str(A_S0(1)),' \eta = ',num2str(alpha/(2*cp_S0(1)))],...
           ['A= ',num2str(A_S0(7)),' \eta = ',num2str(alpha/(2*cp_S0(7)))],...
           ['A= ',num2str(A_S0(10)),' \eta = ',num2str(alpha/(2*cp_S0(10)))],...
           ['A= ',num2str(A_S0(14)),' \eta = ',num2str(alpha/(2*cp_S0(14)))],...
           ['A= ',num2str(A_S0(19)),' \eta = ',num2str(alpha/(2*cp_S0(19)))]);
    title('S0 mode');
    xlim([0 0.2]);ylim([0 0.005]);

% SH0 mode
figure;
    plot(R,Data_polar_SH0_magnitude(1,:),'g');
    hold on;
    plot(R,Data_polar_SH0_magnitude(7,:),'k');
    plot(R,Data_polar_SH0_magnitude(10,:),'r');
    plot(R,Data_polar_SH0_magnitude(14,:),'m');
    plot(R,Data_polar_SH0_magnitude(19,:),'b');

    plot(R,f_SH0(1,:),'g--');
    plot(R,f_SH0(7,:),'k--');
    plot(R,f_SH0(10,:),'r--');
    plot(R,f_SH0(14,:),'m--');
    plot(R,f_SH0(19,:),'b--');

    xlabel('d [m]');
    ylabel('A [m/s]');
    legend([num2str(beta(1)), ' deg'], [num2str(beta(7)), ' deg'], [num2str(beta(10)), ' deg'], ...
           [num2str(beta(14)), ' deg'], [num2str(beta(19)), ' deg'],...
           ['A= ',num2str(A_SH0(1)),' \eta = ',num2str(alpha/(2*cp_SH0(1)))],...
           ['A= ',num2str(A_SH0(7)),' \eta = ',num2str(alpha/(2*cp_SH0(7)))],...
           ['A= ',num2str(A_SH0(10)),' \eta = ',num2str(alpha/(2*cp_SH0(10)))],...
           ['A= ',num2str(A_SH0(14)),' \eta = ',num2str(alpha/(2*cp_SH0(14)))],...
           ['A= ',num2str(A_SH0(19)),' \eta = ',num2str(alpha/(2*cp_SH0(19)))]);
    title('SH0 mode');
    xlim([0 0.2]);ylim([0 0.008]);

% Comparison of amplitudes per mode at 0 deg
figure;
    plot(R,Data_polar_A0_magnitude(1,:),'r');
    hold on;
    plot(R,Data_polar_S0_magnitude(1,:),'b');
    plot(R,Data_polar_SH0_magnitude(1,:),'g');
    legend('A0','S0','SH0');
    title('Amplitude comparison at 0 deg');
    xlabel('d [m]');
    ylabel('A [m/s]');

% Comparison of attenuation per mode at 0 deg
figure;
    plot(R,exp(-eta0.*R)*100,'r--');
    hold on;
    plot(R,exp(-eta0_S0.*R)*100,'b--');
    plot(R,exp(-eta0_SH0.*R)*100,'g--');
    legend('A0','S0','SH0');
    title('Attenuation comparison at 0 deg');
    xlabel('d [m]');
    ylabel('Attenuation [%]');