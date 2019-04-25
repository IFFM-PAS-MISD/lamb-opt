%inputFile 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'weave type';                     fiberType = 'plainWeave';
'stack thickness sequence [m]';   h_p = 0.0980; 
'thickness of the fill [m]';      h_f = 0.048;
'thickness of the warp [m]';      h_w = 0.048;
'width of the fill [m]';          a_f = 0.45;
'thickness of the warp [m]';      a_w = 0.45;
'width of the fill gap [m]';      g_f = 0.3;
'thickness of the warp gap [m]';  g_w = 0.3;   
'volume fraction of fibres';      vol_0 = 0.23;

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
'Young modulus [Pa]';                   e11_f = 275.6;
                                        e22_f = 0.1*e11_f;
'Poisson ratio';                        ni12_f = 0.2;
                                        ni23_f = 0.2;
'density [kg/m3]';                      rho_f = 1900;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'plot weave geometry';                  plot_weave = 'True';
