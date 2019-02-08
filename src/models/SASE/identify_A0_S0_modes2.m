function [FREQ_A0,FREQ_S0,CG_A0,CG_S0] = identify_A0_S0_modes2(FREQ,CG)
% IDENTIFY_A0_S0_MODES   identify A0 and S0 Lamb wave mode dispersion curves
%    extracts only fundamental Lamb wave modes from number_of_modes and
%    angles
% 
% Syntax: [FREQ_A0,FREQ_S0,CG_A0,CG_S0] = identify_A0_S0_modes(FREQ,CG) 
% 
% Inputs: 
%    FREQ - matrix of frequencies, double, dimensions [number_of_modes,number_of_wavenumber_points,number_of_angles], Units: Hz 
%    CG - matrix of group velocities, double, dimensions [number_of_modes,number_of_wavenumber_points,number_of_angles], Units: m/s   
% 
% Outputs: 
%    FREQ_A0 -  matrix of frequencies of A0 mode, double, dimensions [number_of_wavenumber_points,number_of_angles], Units: Hz 
%    FREQ_A0 -  matrix of frequencies of S0 mode, double, dimensions [number_of_wavenumber_points,number_of_angles], Units: Hz 
%    CG_A0 - matrix of group velocities of A0 mode, double, dimensions [number_of_wavenumber_points,number_of_angles], Units: m/s   
%    CG_S0 - matrix of group velocities of S0 mode, double, dimensions [number_of_wavenumber_points,number_of_angles], Units: m/s   
% 
% Example: 
%    [output1,output2] = identify_A0_S0_modes(input1,input2,input3) 
%    [output1,output2] = identify_A0_S0_modes(input1,input2) 
%    [output1] = identify_A0_S0_modes(input1,input2,input3) 
% 
% Other m-files required: none 
% Subfunctions: none 
% MAT-files required: none 
% See also:  
% 

% Author: Pawel Kudela, D.Sc., Ph.D., Eng. 
% Institute of Fluid Flow Machinery Polish Academy of Sciences 
% Mechanics of Intelligent Structures Department 
% email address: pk@imp.gda.pl 
% Website: https://www.imp.gda.pl/en/research-centres/o4/o4z1/people/ 

%---------------------- BEGIN CODE---------------------- 

[~,number_of_wavenumber_points,number_of_angles] = size(FREQ);
FREQ_A0 = zeros(number_of_wavenumber_points,number_of_angles);
FREQ_S0 = zeros(number_of_wavenumber_points,number_of_angles);
CG_A0 = zeros(number_of_wavenumber_points,number_of_angles);
CG_S0 = zeros(number_of_wavenumber_points,number_of_angles);

for j=1:number_of_angles
    f1=FREQ(:,2,j);
    I3=find(f1<50e3); % should find 3 modes existing in range 0-50 kHz, namely A0, S0 and SH0
    FREQ_A0(:,j)=FREQ(I3(1),:,j);
    CG_A0(:,j)=CG(I3(1),:,j);
    if(length(I3)==3)
        FREQ_S0(:,j)=FREQ(I3(3),:,j);
        CG_S0(:,j)=CG(I3(3),:,j);
    else
        FREQ_S0(:,j)=FREQ(I3(2),:,j);
        CG_S0(:,j)=CG(I3(2),:,j);
    end
end

%---------------------- END OF CODE---------------------- 

% ================ [identify_A0_S0_modes.m] ================  
