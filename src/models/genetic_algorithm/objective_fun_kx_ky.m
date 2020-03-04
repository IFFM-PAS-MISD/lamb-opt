function [obj_score] = objective_fun_kx_ky(Data_polar,wavenumber_max,wavenumber,number_of_modes_considered)
% OBJECTIVE_FUN   score for overlap of numerical and experimental dispersion curves
%   kx-ky plane (or rather k-angle plane)
% 
% Syntax: [obj_score] = objective_fun_kx_ky(Data_polar,wavenumber,FREQ,number_of_modes_considered) 
% 
% Inputs: 
%    Data_polar - Experimental dispersion curves, complex double, dimensions [number_of_angles,number_of_wavenumber_points,number_of_frequency_points]  
%    wavenumber - wavenumber matrix computed by SASE mode for frequencies selected from experimental data, double [rad/m]
%    dimensions[number_of_modes,number_of_frequency_points,number_of_angles]
%    FREQ - Numerical frequency matrix, double, 
%           dimensions [number_of_frequency_points,number_of_angles], Units: Hz 
%    number_of_modes_considered - number of modes considered in calculation of the score
% 
% Outputs: 
%    obj_score - Objective function score, double, dimensions [m, n], Units: - 
% 
% Example: 
%    [obj_score] = objective_fun_kx_ky(Data_polar,fmax,FREQ,number_of_modes_considered)
%    [obj_score] = objective_fun_kx_ky(Data_polar,fmax,FREQ) 
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

if nargin == 3
    number_of_modes_considered = 4; % default number of modes considered
end
    
[number_of_angles,number_of_wavenumber_points,number_of_frequency_points] = size(Data_polar);
kvec = linspace(0,wavenumber_max,number_of_wavenumber_points);
obj_score=0;
for j=1:1:number_of_frequency_points 
    %% Function to logical matrix
    % convert numerical dispersion curve into logical matrix (k-beta plane)

    H=zeros(number_of_angles,number_of_wavenumber_points);
    for i=1:number_of_angles % beta
        for k =1:number_of_modes_considered
            [~,I] = min(abs( wavenumber(k,j,i) - kvec )); % mode 1 : number_of_modes_considered (default=4)         
                H(i,I) = 1;
        end
    end
    
    start_idx1 = 1;
    start_idx2 = 3;
    end_idx1 = number_of_angles;
    end_idx2 = number_of_wavenumber_points;
    
    % apply logical matrix to experimental dispersion curve matrix
    dispersion = H(start_idx1:end_idx1,start_idx2:end_idx2).*squeeze(abs(Data_polar(start_idx1:end_idx1,start_idx2:end_idx2,j)));
    
    [I,J] = find(H(start_idx1:end_idx1,start_idx2:end_idx2)==true);

    [m,n]=size(dispersion);
    % sum of values
    for k =1:length(J)
        obj_score = obj_score + ((dispersion(I(k),J(k))))/(m*n);
    end
end

%---------------------- END OF CODE---------------------- 

% ================ [objective_fun_kx_ky.m] ================  
