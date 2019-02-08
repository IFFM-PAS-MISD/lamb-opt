% plot dispersion curves for range of variation

foldername = 'SASE';
modelname = 'SASE1';
data_process_type = 'raw';
data_origin = 'num';
model_output_path = prepare_model_paths(data_process_type,data_origin,foldername,modelname);


test_case = 0;

figure(1)
hold on;
j=1; % beta = 0 deg
for i1 = 1:11
    for i2 = 1:11
        test_case= test_case+1;
        output_name = [model_output_path,filesep,'output',num2str(test_case)];
        load(output_name);
        plot(FREQ_A0(2:end,j)/1e3,wavenumber(2:end,j),'LineWidth',2);% A0
        plot(FREQ_S0(2:end,j)/1e3,wavenumber(2:end,j),'LineWidth',2);% S0
        
    end
end
Hxl=xlabel('Frequency [kHz]');
Hyl=ylabel('Wavenumber [1/m]');
set(Hxl,'FontSize',12);set(Hyl,'FontSize',12);