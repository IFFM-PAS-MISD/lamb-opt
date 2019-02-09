% plot dispersion curves for range of variation
clear all; close all;
foldername = 'SASE';
modelname = 'SASE1';
data_process_type = 'raw';
data_origin = 'num';
model_output_path = prepare_model_paths(data_process_type,data_origin,foldername,modelname);
beta= 0:15:90;

for j=1:length(beta)

    figure(j)
    hold on;
    test_case = 0;
    for i1 = 1:11
        for i2 = 1:11
            test_case= test_case+1;
            output_name = [model_output_path,filesep,'output',num2str(test_case)];
            load(output_name); % FREQ CG wavenumber
            plot(squeeze(FREQ(1,2:end,j))/1e3,wavenumber(2:end,j),'LineWidth',2);% mode 1
            plot(squeeze(FREQ(2,2:end,j))/1e3,wavenumber(2:end,j),'LineWidth',2);% mode 2
            plot(squeeze(FREQ(3,2:end,j))/1e3,wavenumber(2:end,j),'LineWidth',2);% mode 3
            plot(squeeze(FREQ(4,2:end,j))/1e3,wavenumber(2:end,j),'LineWidth',2);% mode 4
            %plot(squeeze(FREQ(5,2:end,j))/1e3,wavenumber(2:end,j),'LineWidth',2);% mode 5
        end
    end
    Hxl=xlabel('f [kHz]');
    Hyl=ylabel('k [1/m]');
    axis([0 600 0 3000]);
    set(Hxl,'FontSize',12);set(Hyl,'FontSize',12);
    title([num2str(beta(j)),' deg']);
%     figfilename = ['reports',filesep,'figures'];
%     print('-dpng', '-r300',figfilename); 
end