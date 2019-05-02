%inputFile 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'weave type';                     fiberType = 'plainWeave';
'lamina thickness [mm]';     h_p = 3/12; % 3 mm total thickness/ 12 layers
'thickness of the fill [mm]';      h_f = h_p/2;
'thickness of the warp [mm]';      h_w = h_p/2;
'width of the fill [mm]';          a_f = 1.92;
'thickness of the warp [mm]';      a_w = 2;
'width of the fill gap [mm]';      g_f = 0;
'thickness of the warp gap [mm]';  g_w = 0;   
'volume fraction of fibres';      vol_0 = 0.5;
%composite material properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matrix(1) - epoxy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'Young modulus [Pa]';                   e11_m = 3.43;
'Poisson ratio';                        ni12_m = 0.35;
'density [kg/m3]';                      rho_m = 1250;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fibres(1) - carbon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'Young modulus [Pa]';                   e11_f = 240;
                                        e22_f = 0.1*e11_f;
'Poisson ratio';                        ni12_f = 0.2;
                                        ni23_f = 0.2;
'density [kg/m3]';                      rho_f = 1900;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'plot weave geometry';                  plot_weave = false;
