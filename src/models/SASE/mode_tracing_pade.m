function [cg_new,om_new] = mode_tracing_pade(cg,om,wavenumber_step)
% MODE_TRACING   Mode tracing using Taylor approximation
%    sorts roots of Lamb waves dispersion curves equations 
% 
% Syntax: [cg_new,om_new] = mode_tracing(cg,om) 
% 
% Inputs: 
%    cg - group velocity dispersion curves, double, dimensions [num_of_modes,number_of_wavenumber_points], Units: m/s 
%    om - angular frequencies, double, dimensions [num_of_modes,number_of_wavenumber_points], Units: [rad/s]
%    wavenumber_step - wavenumber step (equidistant points in wavenumber domain), double, Units: [1/m]
% 
% Outputs: 
%    cg_new - sorted group velocity modes, double, dimensions [num_of_modes,number_of_wavenumber_points], Units: m/s 
%    om_new - sorted angular frequencies, double, double, dimensions [num_of_modes,number_of_wavenumber_points], Units: [rad/s] 
% 
% Example: 
%    [cg_new,om_new] = mode_tracing(cg,om) 
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

[num_of_modes,number_of_wavenumber_points] = size(cg);
ind = zeros(num_of_modes,number_of_wavenumber_points);
cg_new=cg;
om_new=om;
om_max= max(max(om))*10;
for i=1:num_of_modes
    ind(i,number_of_wavenumber_points)=i;
end
for i=1:num_of_modes
    I=i;cc=0;
for k=number_of_wavenumber_points:-1:2
   
    if (cc==0) % first order Taylor approximation
        cc=cc+1;
        cg0=cg(I,k);
        %cg1=cg(I,k-1);
        gamma0=om(I,k); 
        gamma1=cg0;
        om_t=gamma0-gamma1*wavenumber_step; % first order Taylor approximation (backward)
        %cg_t=cg1; 
        %cg_t=cg0; 
        delta_om = abs((om_t-om(:,k-1))/om_t); % deviation of angular frequency in Taylor method
        %delta_cg = abs((cg_t-cg(:,k-1))/cg_t);  % deviation of group velocity in Taylor method
        %delta = delta_om.^2+delta_cg.^2;
        [~,I] = min(delta_om);
        %ind(i,k-1) = I;
        om_new(i,k-1) = om(I,k-1);
        cg_new(i,k-1) = cg(I,k-1);
    elseif (cc==1) % second order Taylor approximation
        cc=cc+1; 
        cg0=cg(I,k+1);
        cg1=cg(I,k);
        gamma0=om(I,k); 
        gamma1=cg0;
        gamma2=1/2*(cg1-cg0)/wavenumber_step;
        om_t=gamma0-gamma1*wavenumber_step-gamma2*wavenumber_step^2; % second order Taylor approximation (backward)
        cg_t=2*cg1-cg0; % linear Taylor approximation
        delta_om = abs((om_t-om(:,k-1))/om_t); % deviation of angular frequency in Taylor method
        delta_cg = abs((cg_t-cg(:,k-1))/cg_t);  % deviation of group velocity in Taylor method
        delta = delta_om.^2+delta_cg.^2;
        [~,I] = min(delta_om);
        %[~,I] = min(delta);
        %ind(i,k-1) = I;
        om_new(i,k-1) = om(I,k-1);
        cg_new(i,k-1) = cg(I,k-1);
        om(I,k) = om_max;
    else
        cc=cc+1; 
        cg0=cg(I,k+2);
        cg1=cg(I,k+1);
        cg2=cg(I,k);
    
        gamma0=om(I,k);
        gamma1=cg2;
        gamma2=1/2*(cg2-cg1)/wavenumber_step;
        gamma3=1/6*(cg2-2*cg1+cg0)/(wavenumber_step^2);
    
        beta2=(gamma2^2-gamma3*gamma1)/(gamma1^2-gamma2*gamma0);
        beta1=-(gamma2+gamma0*beta2)/gamma1;
        alpha0=gamma0;
        alpha1=gamma1+gamma0*beta1;
    
        %om_p=(alpha0 + alpha1 * wavenumber_step)/(1 + beta1*wavenumber_step + beta2*wavenumber_step^2); %Pade expansion of order [1/2](forward)
        om_p=(alpha0 - alpha1 * wavenumber_step)/(1 -( beta1*wavenumber_step + beta2*wavenumber_step^2)); % Pade expansion of order [1/2](backward)
        %om_t=gamma0-gamma1*wavenumber_step-gamma2*wavenumber_step^2-gamma3*wavenumber_step^3; % Third order Taylor approximation (backward)
    
        %cg_t=2*cg2-cg1; % linear Taylor approximation
        
        %delta_om = abs((om_t-om(:,k-1))/om_t); % deviation of angular frequency in Taylor method
        delta_om = abs((om_p-om(:,k-1))/om_p); % deviation of angular frequency in Pade method
        %delta_cg = abs((cg_t-cg(:,k-1))/cg_t);  % deviation of group velocity in Taylor method
        %delta = delta_om.^2+delta_cg.^2;
        %[a,I] = min(delta);
        [~,I] = min(delta_om);
        
        %ind(i,k-1) = I;
        om_new(i,k-1) = om(I,k-1);
        cg_new(i,k-1) = cg(I,k-1);
        om(I,k) = om_max;
        if(I~=i)
            cc=0;
           
        end
    end

   
    
    
end
end 

%---------------------- END OF CODE---------------------- 

% ================ [mode_tracing.m] ================  
