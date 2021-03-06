%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%                        input.m                      %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define material properties here
% rho in kg/m^3
% C in GPa
% And input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Inputs: layup, material properties
% layup angles are as given in the paper: +Z = 0 deg. and then
% counter-clockwise 
% For material properties transformation, +Y = 0 deg. and then
% counter-clockwise -- for stacking direction +X.
% Hence, layup angles have to be changed. rot_angles = layup - 90 deg.

%% Input

% Layup
% layup = [0 0 90 90 90 90 90 90 90 90 0 0];
%layup = [0,0,90,90,0,0];
%layup = [90,90,90,90,90,90];
%layup = [0 90 0 90 0];
% layup = [0 0 0 0 0];
 layup = [0 0 0 0];
%layup = [0 90 0 90 0 90 0 0 90 0 90 0 90 0];
%layup = [0 90 0 90 90 0 90 0];
%layup = [0 90 90 0];
%layup = [0 90 0 90];
%layup = [0 45 90 -45 0];
% layup = [45 0 135 90 90 135 0 45];
% layup = [0 45 -45 90 90 -45 45 0];
%layup = [0];

%layup = [90 90 90 90 90 90 90 90 90 90 90 90]; % S1
%layup = [0 0 0 0 0 0 0 0 0 0 0 0]; % S2
%layup = [0 0 90 90 90 90 90 90 90 90 0 0]; % S3
%layup = [0 90 45 -45 0 90 90 0 -45 45 90 0]; % S4
%layup = [45 135 0 90 90 0 135 45 45 135 0 90 90 0 135 45 45 135 0 90 90 0 135 45];

nlayers = length(layup);
h = 1*ones(nlayers,1) * 1.5e-3/nlayers; % 3 mm thick
% h = 1*ones(nlayers,1) * 3.6*12e-3/12;
%h = 1*ones(nlayers,1) * 1.07*1e-3/8;
%h = 1*ones(nlayers,1) * 2.0*1e-3/12; % [m]
%h = 1*ones(nlayers,1) * 1.56*1e-3/6; % [m]
%h = 1*ones(nlayers,1) * 1.15*1e-3/6; % [m]
%h = 1*ones(nlayers,1) * 3*1e-3/14; % [m]
%h = 1*ones(nlayers,1) * 3.5*1e-3/8; % [m]
%h = 1*ones(nlayers,1) * 1.5*1e-3/4; % [m]
%h = 1*ones(nlayers,1) * 3*1e-3; % [m]
%h = 1*ones(nlayers,1) * 3*1e-3/24; % [m]
%h = 1*ones(nlayers,1) * 1.5*1e-3; % [m]
rot_angles = layup - 90; % set x axis horizontal
    
% Stacking direction
stack_dir = 1;


%% Material Properties

% 1 = Taupin      2 = Neau      3 = Mircea    4 = aluminium 5 = after
% Sierakowski 6=Pan, 7=Short fibers vol 40% 8=short fibers vol 40% 9= short
% fibers vol 20% 10= short fibers vol 10%
props = 5;

%% In Taupin et. al.
switch props
    case 1 
    rho = 1494;
    
    C0r = zeros(6,6);
    C0i = zeros(6,6);
    
    C0r(1,1) = 145.5;
    C0r(2,2) = 11.5;
    C0r(3,3) = 11.5;
    C0r(1,2) = 5.0;
    C0r(1,3) = 5.0;
    C0r(2,3) = 5.258;
    C0r(4,4) = 3.5;
    C0r(5,5) = 5.2;
    C0r(6,6) = 5.2;
    
    C0r(2,1) = C0r(1,2);
    C0r(3,1) = C0r(1,3);
    C0r(3,2) = C0r(2,3);
    
    C0i(1,1) = 10;
    C0i(2,2) = 0.8;
    C0i(3,3) = 0.9;
    C0i(1,2) = 0.7;
    C0i(1,3) = 0.7;
    C0i(2,3) = 0.4;
    C0i(4,4) = 0.2;
    C0i(5,5) = 0.4;
    C0i(6,6) = 0.4;
    
    C0i(2,1) = C0i(1,2);
    C0i(3,1) = C0i(1,3);
    C0i(3,2) = C0i(2,3);
    
    C0r = C0r*1e9;
    C0i = C0i*1e9;
    
    case 2
    % In Neau thesis
    
    rho = 1494;
    
    C0r = zeros(6,6);
    C0i = zeros(6,6);
    
    C0r(1,1) = 86.6;
    C0r(2,2) = 13.5;
    C0r(3,3) = 14.0;
    C0r(1,2) = 9.0;
    C0r(1,3) = 6.4;
    C0r(2,3) = 6.8;
    C0r(4,4) = 2.72;
    C0r(5,5) = 4.06;
    C0r(6,6) = 4.70;
    
    C0r(2,1) = C0r(1,2);
    C0r(3,1) = C0r(1,3);
    C0r(3,2) = C0r(2,3);
    
    C0i(1,1) = 7.50;
    C0i(2,2) = 0.60;
    C0i(3,3) = 0.28;
    C0i(1,2) = 0.30;
    C0i(1,3) = 0.60;
    C0i(2,3) = 0.25;
    C0i(4,4) = 0.10;
    C0i(5,5) = 0.12;
    C0i(6,6) = 0.28;
    
    C0i(2,1) = C0i(1,2);
    C0i(3,1) = C0i(1,3);
    C0i(3,2) = C0i(2,3);
    
    C0r = C0r*1e9;
    C0i = C0i*1e9;
    
    case 3
    
    % In Mircea et. al
    
    rho = 1540;
    
    C0r = zeros(6,6);
    C0i = zeros(6,6);
    
    C0r(1,1) = 125;
    C0r(2,2) = 13.9;
    C0r(3,3) = 14.5;
    C0r(1,2) = 6.3;
    C0r(1,3) = 5.4;
    C0r(2,3) = 7.1;
    C0r(4,4) = 3.7;
    C0r(5,5) = 5.4;
    C0r(6,6) = 5.4;
    
    C0r(2,1) = C0r(1,2);
    C0r(3,1) = C0r(1,3);
    C0r(3,2) = C0r(2,3);
    
    C0i(1,1) = 3;
    C0i(2,2) = 0.6;
    C0i(3,3) = 0.6;
    C0i(1,2) = 0.9;
    C0i(1,3) = 0.4;
    C0i(2,3) = 0.23;
    C0i(4,4) = 0.12;
    C0i(5,5) = 0.3;
    C0i(6,6) = 0.5;
    
    C0i(2,1) = C0i(1,2);
    C0i(3,1) = C0i(1,3);
    C0i(3,2) = C0i(2,3);
    
    C0r = C0r*1e9;
    C0i = C0i*1e9;
    
    case 5
    % calculation according to Sierakowski and standard carbon/epoxy
    rho = 1627;
    
    C0r = zeros(6,6);
    C0i = zeros(6,6);
    
    C0r(1,1) = 163.9;
    C0r(2,2) = 14.14;
    C0r(3,3) = 14.14;
    C0r(1,2) = 5.0;
    C0r(1,3) = 5.0;
    C0r(2,3) = 4.92;
    C0r(4,4) = 4.61;
    C0r(5,5) = 4.6;
    C0r(6,6) = 4.6;
    
    C0r(2,1) = C0r(1,2);
    C0r(3,1) = C0r(1,3);
    C0r(3,2) = C0r(2,3);
    
    C0i(1,1) = 10;
    C0i(2,2) = 0.8;
    C0i(3,3) = 0.9;
    C0i(1,2) = 0.7;
    C0i(1,3) = 0.7;
    C0i(2,3) = 0.4;
    C0i(4,4) = 0.2;
    C0i(5,5) = 0.4;
    C0i(6,6) = 0.4;
    
    C0i(2,1) = C0i(1,2);
    C0i(3,1) = C0i(1,3);
    C0i(3,2) = C0i(2,3);
    
    C0r = C0r*1e9;
    C0i = C0i*1e9;
    
    case 4
    
    % aluminium isotropic case
    
    rho = 2700;
    %rho=rho*0.9983; % 50 Celsius
    C0r = zeros(6,6);
    C0i = zeros(6,6);
    
    C0r(1,1) = 105.2;
    C0r(2,2) = 105.2;
    C0r(3,3) = 105.2;
    C0r(1,2) = 51.8;
    C0r(1,3) = 51.8;
    C0r(2,3) = 51.8;
    C0r(4,4) = 26.7;
    C0r(5,5) = 26.7;
    C0r(6,6) = 26.7;
    
    C0r(2,1) = C0r(1,2);
    C0r(3,1) = C0r(1,3);
    C0r(3,2) = C0r(2,3);
    
    C0i(1,1) = 0;
    C0i(2,2) = 0;
    C0i(3,3) = 0;
    C0i(1,2) = 0;
    C0i(1,3) = 0;
    C0i(2,3) = 0;
    C0i(4,4) = 0;
    C0i(5,5) = 0;
    C0i(6,6) = 0;
    
    C0i(2,1) = C0i(1,2);
    C0i(3,1) = C0i(1,3);
    C0i(3,2) = C0i(2,3);
    
    %C0r = C0r*0.99; % 50 Celsius
    %C0i = C0i*0.99; % 50 Celsius
    
    C0r = C0r*1e9;
    C0i = C0i*1e9;
    
    case 6 % short fibres ROM by Pan
    rho = 1733;
    
    C0r = zeros(6,6);
    C0i = zeros(6,6);
    
    C0r(1,1) = 2.0858864e+01;
    C0r(2,2) = 2.0858864e+01;
    C0r(3,3) = 2.0858864e+01;
    C0r(1,2) = 1.0494701e+01;
    C0r(1,3) = 1.0494701e+01;
    C0r(2,3) = 1.0494701e+01;
    C0r(4,4) = 5.1820819e+00;
    C0r(5,5) = 5.1820819e+00;
    C0r(6,6) = 5.1820819e+00;
    
    C0r(2,1) = C0r(1,2);
    C0r(3,1) = C0r(1,3);
    C0r(3,2) = C0r(2,3);
    
    C0i(1,1) = 0.1;
    C0i(2,2) = 0.1;
    C0i(3,3) = 0.1;
    C0i(1,2) = 0.1;
    C0i(1,3) = 0.1;
    C0i(2,3) = 0.1;
    C0i(4,4) = 0.1;
    C0i(5,5) = 0.1;
    C0i(6,6) = 0.1;
    
    C0i(2,1) = C0i(1,2);
    C0i(3,1) = C0i(1,3);
    C0i(3,2) = C0i(2,3);
    
    C0r = C0r*1e9;
    C0i = C0i*1e9;
    
    case 7 % my averaged model
    rho = 1733;
    
    C0r = zeros(6,6);
    C0i = zeros(6,6);
    
    C0r(1,1) = 2.3261880e+001;
    C0r(2,2) = 2.3296343e+001;
    C0r(3,3) = 1.4009152e+001;
    C0r(1,2) = 8.4323538e+000;
    C0r(1,3) = 5.2709958e+000;
    C0r(1,4) = -4.4951897e-003;
    C0r(2,3) = 5.2712388e+000;
    C0r(2,4) = -6.8088053e-003;
    C0r(3,4) = -7.0911699e-005;
    C0r(4,4) = 4.3840580e+000;%7.4342357e+000;
    C0r(5,5) = 4.3840580e+000;
    C0r(5,6) = 1.8688588e-005;
    C0r(6,6) = 7.4342357e+000;%4.3839740e+000;
    
    C0r(2,1) = C0r(1,2);
    C0r(3,1) = C0r(1,3);
    C0r(3,2) = C0r(2,3);
    C0r(4,1) = C0r(1,4);
    C0r(4,2) = C0r(2,4);
    C0r(4,3) = C0r(3,4);
    C0r(6,5) = C0r(5,6);
    
    C0i(1,1) = 0.1;
    C0i(2,2) = 0.1;
    C0i(3,3) = 0.1;
    C0i(1,2) = 0.1;
    C0i(1,3) = 0.1;
    C0i(1,4) = 0.1;
    C0i(2,3) = 0.1;
    C0i(2,4) = 0.1;
    C0i(3,4) = 0.1;
    C0i(4,4) = 0.1;
    C0i(5,5) = 0.1;
    C0i(5,6) = 0.1;
    C0i(6,6) = 0.1;
    
    C0i(2,1) = C0i(1,2);
    C0i(3,1) = C0i(1,3);
    C0i(3,2) = C0i(2,3);
    C0i(4,1) = C0i(1,4);
    C0i(4,2) = C0i(2,4);
    C0i(4,3) = C0i(3,4);
    C0i(6,5) = C0i(5,6);
    
    C0r = C0r*1e9;
    C0i = C0i*1e9;
    
    case 8 % my averaged model v2
    rho = 1733;
    
    C0r = zeros(6,6);
    C0i = zeros(6,6);
    
    C0r(1,1) = 2.3261880e+001;
    C0r(2,2) = 2.3296343e+001;
    C0r(3,3) = 1.4009152e+001;
    C0r(1,2) = 8.4323538e+000;
    C0r(1,3) = 5.2709958e+000;
    C0r(1,4) = 0;
    C0r(2,3) = 5.2712388e+000;
    C0r(2,4) = 0;
    C0r(3,4) = 0;
    C0r(4,4) = 4.3840580e+000;
    C0r(5,5) = 4.3840580e+000;
    C0r(5,6) = 1.8688588e-005;
    C0r(6,6) = 7.4342357e+000;
    
    C0r(2,1) = C0r(1,2);
    C0r(3,1) = C0r(1,3);
    C0r(3,2) = C0r(2,3);
    C0r(4,1) = C0r(1,4);
    C0r(4,2) = C0r(2,4);
    C0r(4,3) = C0r(3,4);
    C0r(6,5) = C0r(5,6);
    
    C0i(1,1) = 0.1;
    C0i(2,2) = 0.1;
    C0i(3,3) = 0.1;
    C0i(1,2) = 0.1;
    C0i(1,3) = 0.1;
    C0i(1,4) = 0;
    C0i(2,3) = 0.1;
    C0i(2,4) = 0;
    C0i(3,4) = 0;
    C0i(4,4) = 0.1;
    C0i(5,5) = 0.1;
    C0i(5,6) = 0.1;
    C0i(6,6) = 0.1;
    
    C0i(2,1) = C0i(1,2);
    C0i(3,1) = C0i(1,3);
    C0i(3,2) = C0i(2,3);
    C0i(4,1) = C0i(1,4);
    C0i(4,2) = C0i(2,4);
    C0i(4,3) = C0i(3,4);
    C0i(6,5) = C0i(5,6);
    
    C0r = C0r*1e9;
    C0i = C0i*1e9;
    
    case 9 % my averaged model v2 vol20%
    rho = 1733;
    
    C0r = zeros(6,6);
    C0i = zeros(6,6);
    
    C0r(1,1) = 1.4540676e+001;
    C0r(2,2) = 1.4625409e+001;
    C0r(3,3) = 9.8547175e+000;
    C0r(1,2) = 5.7816322e+000;
    C0r(1,3) = 4.1481537e+000;
    C0r(1,4) = 0;
    C0r(2,3) = 4.1487257e+000;
    C0r(2,4) = 0;
    C0r(3,4) = 0;
    C0r(4,4) = 2.8608669e+000;
    C0r(5,5) = 2.8610477e+000;
    C0r(5,6) = 6.5375956e-005;
    C0r(6,6) = 4.4379015e+000;
    
    C0r(2,1) = C0r(1,2);
    C0r(3,1) = C0r(1,3);
    C0r(3,2) = C0r(2,3);
    C0r(4,1) = C0r(1,4);
    C0r(4,2) = C0r(2,4);
    C0r(4,3) = C0r(3,4);
    C0r(6,5) = C0r(5,6);
    
    C0i(1,1) = 0.1;
    C0i(2,2) = 0.1;
    C0i(3,3) = 0.1;
    C0i(1,2) = 0.1;
    C0i(1,3) = 0.1;
    C0i(1,4) = 0;
    C0i(2,3) = 0.1;
    C0i(2,4) = 0;
    C0i(3,4) = 0;
    C0i(4,4) = 0.1;
    C0i(5,5) = 0.1;
    C0i(5,6) = 0.1;
    C0i(6,6) = 0.1;
    
    C0i(2,1) = C0i(1,2);
    C0i(3,1) = C0i(1,3);
    C0i(3,2) = C0i(2,3);
    C0i(4,1) = C0i(1,4);
    C0i(4,2) = C0i(2,4);
    C0i(4,3) = C0i(3,4);
    C0i(6,5) = C0i(5,6);
    
    C0r = C0r*1e9;
    C0i = C0i*1e9;
    
    case 10 % my averaged model v2 vol10%
    rho = 1733;
    
    C0r = zeros(6,6);
    C0i = zeros(6,6);
    
    C0r(1,1) = 9.8436725e+000;
    C0r(2,2) = 9.7539428e+000;
    C0r(3,3) = 7.5637467e+000;
    C0r(1,2) = 4.2848850e+000;
    C0r(1,3) = 3.5295627e+000;
    C0r(1,4) = 0;
    C0r(2,3) = 3.5289642e+000;
    C0r(2,4) = 0;
    C0r(3,4) = 0;
    C0r(4,4) = 2.0208502e+000;
    C0r(5,5) = 2.0210344e+000;
    C0r(5,6) = -2.4745421e-005;
    C0r(6,6) = 2.7502190e+000;
    
    C0r(2,1) = C0r(1,2);
    C0r(3,1) = C0r(1,3);
    C0r(3,2) = C0r(2,3);
    C0r(4,1) = C0r(1,4);
    C0r(4,2) = C0r(2,4);
    C0r(4,3) = C0r(3,4);
    C0r(6,5) = C0r(5,6);
    
    C0i(1,1) = 0.1;
    C0i(2,2) = 0.1;
    C0i(3,3) = 0.1;
    C0i(1,2) = 0.1;
    C0i(1,3) = 0.1;
    C0i(1,4) = 0;
    C0i(2,3) = 0.1;
    C0i(2,4) = 0;
    C0i(3,4) = 0;
    C0i(4,4) = 0.1;
    C0i(5,5) = 0.1;
    C0i(5,6) = 0.1;
    C0i(6,6) = 0.1;
    
    C0i(2,1) = C0i(1,2);
    C0i(3,1) = C0i(1,3);
    C0i(3,2) = C0i(2,3);
    C0i(4,1) = C0i(1,4);
    C0i(4,2) = C0i(2,4);
    C0i(4,3) = C0i(3,4);
    C0i(6,5) = C0i(5,6);
    
    C0r = C0r*1e9;
    C0i = C0i*1e9;
    
    case 11 
    rho = 1494;
    
    C0r = zeros(6,6);
    C0i = zeros(6,6);
    
    C0r(1,1) = 165;
    C0r(2,2) = 14.5;
    C0r(3,3) = 14.5;
    C0r(1,2) = 6.8;
    C0r(1,3) = 6.8;
    C0r(2,3) = 6.49;
    C0r(4,4) = 3.9;
    C0r(5,5) = 5.2;
    C0r(6,6) = 5.2;
    
    C0r(2,1) = C0r(1,2);
    C0r(3,1) = C0r(1,3);
    C0r(3,2) = C0r(2,3);
    
    C0i(1,1) = 10;
    C0i(2,2) = 0.8;
    C0i(3,3) = 0.9;
    C0i(1,2) = 0.7;
    C0i(1,3) = 0.7;
    C0i(2,3) = 0.4;
    C0i(4,4) = 0.2;
    C0i(5,5) = 0.4;
    C0i(6,6) = 0.4;
    
    C0i(2,1) = C0i(1,2);
    C0i(3,1) = C0i(1,3);
    C0i(3,2) = C0i(2,3);
    
    C0r = C0r*1e9;
    C0i = C0i*1e9;
end


%%
% save MatProps.mat C0r C0i

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%