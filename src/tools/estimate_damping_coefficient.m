function [alpha_r] = estimate_damping_coefficient(Data,time,alpha)
% ESTIMATE_DAMPING_COEFFICIENT   Estimate damping coefficient based on full wavefield data 
%    Algorithm is based on normalized energy of guided waves 
%    It uses simple grid search over damping coefficients alpha
%    and fits exponent function in the form f(t) = exp(-alpha*t) 
% 
% Syntax: [alpha_r] = estimate_damping_coefficient(Data,time,alpha) 
% 
% Inputs: 
%    Data - Wavefield data, double matrix, dimensions [nx,ny,nft]
%    nx - number of points in x direction
%    ny - number of points in y direction
%    nft - number of time steps
%    time - time vector, dimensions [1,nft], Units: s
%    alpha - optional vector of damping coefficients to fit
% 
% Outputs: 
%    alpha_r - resulting damping coefficient, double, Units: - 
% 
% Example: 
%    [alpha_r] = estimate_damping_coefficient(Data,time,alpha)
%    [alpha_r] = estimate_damping_coefficient(Data,time) 
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

% check arguments
if nargin == 2
    alpha=[1e3:1:1e5]'; % alpha range
end
% calculate energy
%[nx,ny,nft] = size(Data);
%E = zeros(nft,1);
%for k=1:nft
%    for i=1:nx
%        for j=1:ny
%         E(k)=E(k)+Data(i,j,k)^2;
%        end
%    end
%end
Data = Data.^2;
E= squeeze(sum(sum(Data,1),2));
En = E/max(E); % normalization

%% estimate damping
[A,I] = max(En);
Ens= En(I:end);
t=linspace(0,time(length(Ens)),length(Ens));
 
% grid search optimization
f=exp(-alpha*t);
Ens_p = repmat(Ens',[length(alpha),1]);
f2=sum((f-Ens_p).^2,2); 
[~,J] = min(f2);

% sanity check - plot results
figure;
plot(t,f(J,:),'r');
hold on; plot(t,Ens,'k','LineWidth',2);
set(gcf,'color','white');
xlabel('Time [s]');
ylabel('Normalized energy');

alpha_r=alpha(J);


%---------------------- END OF CODE---------------------- 

% ================ [estimate_damping_coefficient.m] ================  
