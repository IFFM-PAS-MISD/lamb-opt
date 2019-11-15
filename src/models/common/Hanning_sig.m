function [t,st] = Hanning_sig(dt,nft,fm,fc,t_1,w)
% HANNING_SIG   One line description of what the function or script performs (H1 line) 
%    optional: more details about the function than in the H1 line 
%    optional: more details about the function than in the H1 line 
%    optional: more details about the function than in the H1 line 
% 
% Syntax: [output1,output2] = Hanning_sig(input1,input2,input3) 
% 
% Inputs: 
%    input1 - Description, string, dimensions [m, n], Units: ms 
%    input2 - Description, logical, dimensions [m, n], Units: m 
%    input3 - Description, double, dimensions [m, n], Units: N 
% 
% Outputs: 
%    output1 - Description, integer, dimensions [m, n], Units: - 
%    output2 - Description, double, dimensions [m, n], Units: m/s^2 
% 
% Example: 
%    [output1,output2] = Hanning_sig(input1,input2,input3) 
%    [output1,output2] = Hanning_sig(input1,input2) 
%    [output1] = Hanning_sig(input1,input2,input3) 
% 
% Other m-files required: none 
% MAT-files required: none 
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2 
% 

% Author: Pawel Kudela, D.Sc., Ph.D., Eng. 
% Institute of Fluid Flow Machinery Polish Academy of Sciences 
% Mechanics of Intelligent Structures Department 
% email address: pk@imp.gda.pl 
% Website: https://www.imp.gda.pl/en/research-centres/o4/o4z1/people/ 

%---------------------- BEGIN CODE---------------------- 

t_t=1/fm;   % total duration time of the excitation [s]
t_2=t_1+t_t; % excitation termiation time [s]
t=[0:nft-1]*dt; % time vector

st=zeros(nft,1);

for nf=nft/2:nft/2+w+1
    st(nf)=0.5*(1-cos(2*pi*fm*(t(nf)-t_2)))*sin(2*pi*fc*(t(nf)-t_1));
end

%---------------------- END OF CODE---------------------- 

% ================ [Hanning_sig.m] ================  
