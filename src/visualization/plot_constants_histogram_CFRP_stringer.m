clear all; close all;
% load projectroot path
load project_paths projectroot src_path;
overwrite = true; % allow overwriting existing results if true

% figure parameters
% size 12cm by 8cm (1-column text)
%fig_width = 12; fig_height = 8; 
% size 7cm by 5cm (2-column text)
fig_width = 7; fig_height = 5; 
% create path to the numerical model data folder
modelfolder = 'genetic_algorithm';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% direct method   k-f -90deg
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 modelname1 = 'ga_stringer_C_tensor_known_mass';
% create path to the numerical raw data folder
model_input_path = fullfile( projectroot, 'data','raw','num',modelfolder,[modelname1,'_out'], filesep );


C11_=zeros(100,1);
C12_=zeros(100,1);
C13_=zeros(100,1);
C22_=zeros(100,1);
C23_=zeros(100,1);
C33_=zeros(100,1);
C44_=zeros(100,1);
C55_=zeros(100,1);
C66_=zeros(100,1);

ObjVal_=zeros(100,1);
c=0;
range=[1:10];

for test_case=range
    c=c+1;
    output_name = [model_input_path,filesep,num2str(test_case),'output'];
    load(output_name); % 'C11','C12','C13','C22','C23','C33','C44','C55','C66','ObjVal'

    C11_(test_case) = C11;
    C12_(test_case) = C12;
    C13_(test_case) = C13;
    C22_(test_case) = C22;
    C23_(test_case) = C23;
    C33_(test_case) = C33;
    C44_(test_case) = C44;
    C55_(test_case) = C55;
    C66_(test_case) = C66;
    ObjVal_(test_case)=ObjVal;
end

C11mean = mean(C11_(range));
C11std = std(C11_(range));
C12mean = mean(C12_(range));
C12std = std(C12_(range));
C13mean = mean(C13_(range));
C13std = std(C13_(range));
C22mean = mean(C22_(range));
C22std = std(C22_(range));
C23mean = mean(C23_(range));
C23std = std(C23_(range));
C33mean = mean(C33_(range));
C33std = std(C33_(range));
C44mean = mean(C44_(range));
C44std = std(C44_(range));
C55mean = mean(C55_(range));
C55std = std(C55_(range));
C66mean = mean(C66_(range));
C66std = std(C66_(range));

[Fd,I]=min(ObjVal_(range));
Fdmean=mean(ObjVal_(range));
Fdstd=std(ObjVal_(range));

C11best=C11_(range(I));
C12best=C12_(range(I));
C13best=C13_(range(I));
C22best=C22_(range(I));
C23best=C23_(range(I));
C33best=C33_(range(I));
C44best=C44_(range(I));
C55best=C55_(range(I));
C66best=C66_(range(I));  

  Tab_direct = [C11best, C11mean, C11std;
                          C12best, C12mean, C12std;
                          C13best, C13mean, C13std;
                          C22best, C22mean, C22std;
                          C23best, C23mean, C23std;
                          C33best, C33mean, C33std;
                          C44best, C44mean, C44std;
                          C55best, C55mean, C55std;
                          C66best, C66mean, C66std]/1e9;
 
 Cdbest=cellstr(num2str([Tab_direct(:,1);Fd],'%.2f'));
 Cdmean=cellstr(num2str([Tab_direct(:,2);Fdmean],'%.2f'));
 Cdstd=cellstr(num2str([Tab_direct(:,3);Fdstd],'%.2f'));
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% direct method  k-f 0 deg
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 modelname2 = 'ga_stringer_C_tensor_known_mass1';
% create path to the numerical raw data folder
model_input_path = fullfile( projectroot, 'data','raw','num',modelfolder,[modelname2,'_out'], filesep );


C11_=zeros(100,1);
C12_=zeros(100,1);
C13_=zeros(100,1);
C22_=zeros(100,1);
C23_=zeros(100,1);
C33_=zeros(100,1);
C44_=zeros(100,1);
C55_=zeros(100,1);
C66_=zeros(100,1);

ObjVal_=zeros(100,1);
c=0;
range=[1];

for test_case=range
    c=c+1;
    output_name = [model_input_path,filesep,num2str(test_case),'output'];
    load(output_name); % 'C11','C12','C13','C22','C23','C33','C44','C55','C66','ObjVal'

    C11_(test_case) = C11;
    C12_(test_case) = C12;
    C13_(test_case) = C13;
    C22_(test_case) = C22;
    C23_(test_case) = C23;
    C33_(test_case) = C33;
    C44_(test_case) = C44;
    C55_(test_case) = C55;
    C66_(test_case) = C66;
    ObjVal_(test_case)=ObjVal;
end

C11mean = mean(C11_(range));
C11std = std(C11_(range));
C12mean = mean(C12_(range));
C12std = std(C12_(range));
C13mean = mean(C13_(range));
C13std = std(C13_(range));
C22mean = mean(C22_(range));
C22std = std(C22_(range));
C23mean = mean(C23_(range));
C23std = std(C23_(range));
C33mean = mean(C33_(range));
C33std = std(C33_(range));
C44mean = mean(C44_(range));
C44std = std(C44_(range));
C55mean = mean(C55_(range));
C55std = std(C55_(range));
C66mean = mean(C66_(range));
C66std = std(C66_(range));

[Fd,I]=min(ObjVal_(range));
Fdmean=mean(ObjVal_(range));
Fdstd=std(ObjVal_(range));

C11best=C11_(range(I));
C12best=C12_(range(I));
C13best=C13_(range(I));
C22best=C22_(range(I));
C23best=C23_(range(I));
C33best=C33_(range(I));
C44best=C44_(range(I));
C55best=C55_(range(I));
C66best=C66_(range(I));  

  Tab_direct = [C11best, C11mean, C11std;
                          C12best, C12mean, C12std;
                          C13best, C13mean, C13std;
                          C22best, C22mean, C22std;
                          C23best, C23mean, C23std;
                          C33best, C33mean, C33std;
                          C44best, C44mean, C44std;
                          C55best, C55mean, C55std;
                          C66best, C66mean, C66std]/1e9;
 
 Cebest=cellstr(num2str([Tab_direct(:,1);Fd],'%.2f'));
 Cemean=cellstr(num2str([Tab_direct(:,2);Fdmean],'%.2f'));
 Cestd=cellstr(num2str([Tab_direct(:,3);Fdstd],'%.2f'));
 %% plot histogram
figure;
hist(C11_(range),8);
Models=[modelname1,'_',modelname2];
 % prepare csv table for latex document
 paper_path=[projectroot,'reports',filesep,'journal_papers',filesep,'Composite_Structures_unidirectional',filesep];
 Constants = {'$C_{11}$';'$C_{12}$';'$C_{13}$';'$C_{22}$';'$C_{23}$';'$C_{33}$';'$C_{44}$';'$C_{55}$';'$C_{66}$';'$F$'};
 T = table(Cdbest,Cdmean,Cdstd,Cebest,Cemean,Cestd,'RowNames',Constants);
 writetable(T,[paper_path,Models,'.csv'],'WriteRowNames',true);
 