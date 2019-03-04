% plot objective function score (bar plot)

clear all; close all;

% load projectroot path
load project_paths projectroot src_path;
pallete = 'high-contrast';
run(fullfile(projectroot,'src','tools','colours.m'));
overwrite = false; % allow overwriting existing results if true

% figure parameters
% size 12cm by 8cm (1-column text)
fig_width = 12; fig_height = 8; 
% size 7cm by 5cm (2-column text)
%fig_width = 7; fig_height = 5; 
% create path to the numerical model data folder
modelfolder = 'SASE';
modelname = 'SASE1';
radians = false;
% create output path
output_path = prepare_figure_paths(modelfolder,modelname);

% create path to the numerical raw data folder
model_input_path = fullfile( projectroot, 'data','interim','num',modelfolder,[modelname,'_out'], filesep );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 MODE CONSIDERED (1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([model_input_path,'objective_fun_score_nmodes_1']); %objective_fun_score
load([model_input_path,'objective_fun_simul_pzt_score_nmodes_1']); % objective_fun_pzt_score
[~,I1] = max(objective_fun_score);
[~,I2] = max(objective_fun_pzt_score);
y = [objective_fun_score/max(objective_fun_score),objective_fun_pzt_score/max(objective_fun_pzt_score)];

figure(1);box on;
b=bar(y,1,'FaceColor','flat');

for j = 1:size(y,1)
    b(1).CData(j,:) = YELLOW/255; % first bar group color
end
for j = 1:size(y,1)
    b(2).CData(j,:) = RED/255; % second bar group color
end
b(1).CData(I1,:) = BLUE/255; % best fit for laser data
b(2).CData(I2,:) = BLACK/255; % best fit for pzt data

hold on;
bar(NaN,'FaceColor',BLUE/255,'EdgeColor','none');
bar(NaN,'FaceColor',BLACK/255,'EdgeColor','none');
l = cell(1,4);
l{1}='laser'; l{2}='pzt'; l{3}='laser opt'; l{4}='pzt opt';
legend(l,'Location','SouthWest');

set(gca,'FontName','Times');
set(gcf,'Color','w');
set(gca, 'Layer', 'Top');
set(gca,'Fontsize',10,'linewidth',1);

fig = gcf;
title(['Score for mode 1, laser opt = ',num2str(I1),', pzt opt = ',num2str(I2)],'Fontsize',12,'interpreter','latex');
set(fig, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); % size 12cm by 8cm (1-column text)
% remove unnecessary white space
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
fig.PaperPositionMode   = 'auto';
figfilename = 'objective_function_score_modes_1_';
print([output_path,figfilename],'-dpng', '-r600'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2 MODES CONSIDERED (1,2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([model_input_path,'objective_fun_score_nmodes_2']); %objective_fun_score
load([model_input_path,'objective_fun_simul_pzt_score_nmodes_2']); % objective_fun_pzt_score
[~,I1] = max(objective_fun_score);
[~,I2] = max(objective_fun_pzt_score);
y = [objective_fun_score/max(objective_fun_score),objective_fun_pzt_score/max(objective_fun_pzt_score)];

figure(2);box on;
b=bar(y,1,'FaceColor','flat');

for j = 1:size(y,1)
    b(1).CData(j,:) = YELLOW/255; % first bar group color
end
for j = 1:size(y,1)
    b(2).CData(j,:) = RED/255; % second bar group color
end
b(1).CData(I1,:) = BLUE/255; % best fit for laser data
b(2).CData(I2,:) = BLACK/255; % best fit for pzt data

hold on;
bar(NaN,'FaceColor',BLUE/255,'EdgeColor','none');
bar(NaN,'FaceColor',BLACK/255,'EdgeColor','none');
l = cell(1,4);
l{1}='laser'; l{2}='pzt'; l{3}='laser opt'; l{4}='pzt opt';
legend(l,'Location','SouthWest');

set(gca,'FontName','Times');
set(gcf,'Color','w');
set(gca, 'Layer', 'Top');
set(gca,'Fontsize',10,'linewidth',1);
fig = gcf;
title(['Score for modes 1-2, laser opt = ',num2str(I1),', pzt opt = ',num2str(I2)],'Fontsize',12,'interpreter','latex');
set(fig, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); % size 12cm by 8cm (1-column text)
% remove unnecessary white space
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
fig.PaperPositionMode   = 'auto';
figfilename = 'objective_function_score_modes_2_';
print([output_path,figfilename],'-dpng', '-r600'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3 MODES CONSIDERED (1,2,3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([model_input_path,'objective_fun_score_nmodes_3']); %objective_fun_score
load([model_input_path,'objective_fun_simul_pzt_score_nmodes_3']); % objective_fun_pzt_score
[~,I1] = max(objective_fun_score);
[~,I2] = max(objective_fun_pzt_score);
y = [objective_fun_score/max(objective_fun_score),objective_fun_pzt_score/max(objective_fun_pzt_score)];

figure(3);box on;
b=bar(y,1,'FaceColor','flat');

for j = 1:size(y,1)
    b(1).CData(j,:) = YELLOW/255; % first bar group color
end
for j = 1:size(y,1)
    b(2).CData(j,:) = RED/255; % second bar group color
end
b(1).CData(I1,:) = BLUE/255; % best fit for laser data
b(2).CData(I2,:) = BLACK/255; % best fit for pzt data

hold on;
bar(NaN,'FaceColor',BLUE/255,'EdgeColor','none');
bar(NaN,'FaceColor',BLACK/255,'EdgeColor','none');
l = cell(1,4);
l{1}='laser'; l{2}='pzt'; l{3}='laser opt'; l{4}='pzt opt';
legend(l,'Location','SouthWest');

set(gca,'FontName','Times');
set(gcf,'Color','w');
set(gca, 'Layer', 'Top');
set(gca,'Fontsize',10,'linewidth',1);
fig = gcf;
title(['Score for modes 1-3, laser opt = ',num2str(I1),', pzt opt = ',num2str(I2)],'Fontsize',12,'interpreter','latex');
set(fig, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); % size 12cm by 8cm (1-column text)
% remove unnecessary white space
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
fig.PaperPositionMode   = 'auto';
figfilename = 'objective_function_score_modes_3_';
print([output_path,figfilename],'-dpng', '-r600'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4 MODES CONSIDERED (1,2,3,4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([model_input_path,'objective_fun_score_nmodes_4']); %objective_fun_score
load([model_input_path,'objective_fun_simul_pzt_score_nmodes_4']); % objective_fun_pzt_score
[~,I1] = max(objective_fun_score);
[~,I2] = max(objective_fun_pzt_score);
y = [objective_fun_score/max(objective_fun_score),objective_fun_pzt_score/max(objective_fun_pzt_score)];

figure(4);box on;
b=bar(y,1,'FaceColor','flat');

for j = 1:size(y,1)
    b(1).CData(j,:) = YELLOW/255; % first bar group color
end
for j = 1:size(y,1)
    b(2).CData(j,:) = RED/255; % second bar group color
end
b(1).CData(I1,:) = BLUE/255; % best fit for laser data
b(2).CData(I2,:) = BLACK/255; % best fit for pzt data

hold on;
bar(NaN,'FaceColor',BLUE/255,'EdgeColor','none');
bar(NaN,'FaceColor',BLACK/255,'EdgeColor','none');
l = cell(1,4);
l{1}='laser'; l{2}='pzt'; l{3}='laser opt'; l{4}='pzt opt';
legend(l,'Location','SouthWest');

set(gca,'FontName','Times');
set(gcf,'Color','w');
set(gca, 'Layer', 'Top');
set(gca,'Fontsize',10,'linewidth',1);
fig = gcf;
title(['Score for modes 1-4, laser opt = ',num2str(I1),', pzt opt = ',num2str(I2)],'Fontsize',12,'interpreter','latex');
set(fig, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); % size 12cm by 8cm (1-column text)
% remove unnecessary white space
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
fig.PaperPositionMode   = 'auto';
figfilename = 'objective_function_score_modes_4_';
print([output_path,figfilename],'-dpng', '-r600'); 

close all;