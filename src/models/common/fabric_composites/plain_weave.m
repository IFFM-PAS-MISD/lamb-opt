
function [Hx_1,Hx_2,Hx_3,Hx_4, Theta_f, Theta_w, vol_fw] = plain_weave(...
   h, h_f, h_w, a_f, a_w, g_f, g_w, v_f0, nNodes,plot_weave)
% COMPFABRICPROP   Geometry of plain weave 
%    optional: more details about the function than in the H1 line 
%    optional: more details about the function than in the H1 line 
%    optional: more details about the function than in the H1 line 
% 
% Syntax: [output1,output2] = compfabricprop(input1,input2,input3) 
% 
% Inputs: 
%    h - ply thickness, dimensions [1, 1], Units: m 
%    h_f - fill height, dimensions [1, 1], Units: m
%    h_w - warp height, dimensions [1, 1], Units: m
%    a_f - fill width, dimensions [1, 1], Units: m
%    a_w - warp width, dimensions [1, 1], Units: m
%    g_f - gap width along the fill, dimensions [1, 1], Units: m
%    g_w - gap width along the warp, dimensions [1, 1], Units: m
%    vol_0 - overall fiber volume fraction , dimensions [1, 1], Units: -
%    plot_weave - plot weave geometry if true, logical, troue or false, [-]
% 
% Outputs: 
%    Hx_1 - height of matrix 1 in RVE, dimensions [nNodes, nNodes], Units: m
%    Hx_2 - height of fill in RVE, dimensions [nNodes, nNodes], Units: m 
%    Hx_3 - height of warp in RVE, dimensions [nNodes, nNodes], Units: m 
%    Hx_4 - height of matrix 2 in RVE, dimensions [nNodes, nNodes], Units: m 
%    Theta_f - fill undulation angle, dimensions [nNodes, nNodes], Units:
%    rad
%    Theta_w - warp undulation angle, dimensions [nNodes, nNodes], Units:
%    rad
%    vol_fw - fiber volume fraction in fill and warp, dimensions [1, 1], Units: -    
% 
% Example: 
%    [output1,output2] = plain_weave(input1,input2,input3) 
% 
% Other m-files required: 
% Subfunctions: none 
% MAT-files required: none 
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2 
% 

% Author: Piotr Fiborek, D.Sc., Eng. 
% Institute of Fluid Flow Machinery Polish Academy of Sciences 
% Mechanics of Intelligent Structures Department 
% email address: pfiborek@imp.gda.pl 
% Website: https://www.imp.gda.pl/en/research-centres/o4/o4z1/people/ 

%---------------------- BEGIN CODE----------------------
h_m = h-h_f-h_w; % height of matrix

x = linspace(0,(a_w+g_w)/2,nNodes); y = linspace(0,(a_f+g_f)/2,nNodes);
[X,Y] = meshgrid(x,y);


z_yt = h_f/2*cos(pi*a_f/2/(a_f+g_f));
a_yt = pi*a_f/2/(pi-acos(2*z_yt/h_f));
zy_1 = @(y) -h_f/2*cos(pi*y/a_yt);
zy_2 = @(y) h_f/2*cos(pi*y/(a_f+g_f));

 
hy_1 = @(y) (h_f+h_m)/2-zy_2(y);
hy_2 = h_w;
hy_3 = @(y) zy_2(y)-zy_1(y);      % for y=[0,a_f/2]
hy_4 = @(y) (h_f+h_w)/2-zy_1(y); 


z_xt = -h_w/2*cos(pi*a_w/2/(a_w+g_w));
a_xt = pi*a_w/2/acos(2*z_xt/h_w);
zx_1 = @(x,y) h_w/2*cos(pi*x/a_xt)-hy_1(y)+h_m/2;
zx_2 = @(x,y) -h_w/2*cos(pi*x/(a_w+g_w))-hy_1(y)+h_m/2;

hx_1 = @(x,y) (h_w+h_m)/2-zx_1(x,y);
hx_2 = @(x,y) zx_1(x,y) - zx_2(x,y); % for x=[0,a_w/2]
hx_3 = @(x,y) hy_3(y);
hx_4 = @(x,y) zx_2(x,y)-hx_3(x,y)+(h_w+h_m)/2+h_f;

Hx_1 = zeros(size(X));
Hx_1(X<=a_w/2) = hx_1(X(X<=a_w/2),Y(X<=a_w/2));
Hx_1(X>a_w/2) = (h_w+h_m)/2-zx_2(X(X>a_w/2),Y(X>a_w/2));


Hx_2 = hx_2(X,Y);
Hx_2(X > a_w/2 & X <= a_w/2+g_w/2) = 0;

Hx_3 = hx_3(X,Y);
Hx_3(Y > a_f/2 & Y <= a_f/2+g_f/2) = 0;


Hx_4 = zeros(size(X));
Hx_4(Y<=a_f/2) = hx_4(X(Y<=a_f/2),Y(Y<=a_f/2));
Hx_4(Y>a_f/2) = h-(h_w+h_m)/2+zx_2(X(Y>a_f/2),Y(Y>a_f/2));

% variable thickness of the tows cross-section
ef = @(y) abs(zy_2(y)-zy_1(y));              % for y=[0,a_f/2]
ew = @(x) abs(zx_1(x,0)-zx_2(x,0));              % for x=[0,a_w/2]
   
% undulation angle - abs(atan(d/dx zx_2(x,y))), 

thetaf = @(x) abs(atan(pi*h_w/2/(a_w+g_w)*sin(pi*x/(a_w+g_w))));
thetaw = @(y) abs(atan(pi*h_f/2/(a_f+g_f)*sin(pi*y/(a_f+g_f))));

Theta_f = thetaf(X);
Theta_w = thetaw(Y);

% cross-sectional area of the fill and warp tows
Af = quadl(ef,0,a_f/2);
Aw = quadl(ew,0,a_w/2);

% developed lengths of the undulated fill and warp tows inside of the RVE
varf_int = @(x) sqrt(1+(pi*h_w/2/(a_w+g_w)*sin(pi*x/(a_w+g_w))).^2);
Lf = quadl(varf_int,0,(a_w+g_w)/2);

varw_int = @(y) sqrt(1+(pi*h_f/2/(a_f+g_f)*sin(pi*y/(a_f+g_f))).^2);
Lw = quadl(varw_int,0,(a_f+g_f)/2);

% volumes occupied by the fill and warp tows inside of the RVE 

vf = Af*Lf;
vw = Aw*Lw;

% fiber volume fraction inside the fill and warp tows
vol_fw = v_f0*(h*(a_f+g_f)*(a_w+g_w))/4/(vf+vw);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_weave
    x = linspace(0,(a_w+g_w)/2,nNodes); y = linspace(0,(a_f+g_f)/2,nNodes);
    figure(1)
    hold on
        x1 = x(x<=a_w/2); y1=y;
        [X1,Y1] = meshgrid(x1,y1);
        
        surf(X1,Y1,zx_1(X1,Y1))
        surf(X1,Y1,zx_2(X1,Y1))
             
        x2 = x; y2 = y(y <= a_f/2);
        [X2,Y2] = meshgrid(x2,y2);
        surf(X2,Y2,zx_2(X2,Y2))
        surf(X2,Y2,(zx_2(X2,Y2)-hx_3(X2,Y2)))
    shading interp
end
        
       



%---------------------- END OF CODE---------------------- 

% ================ [compfabricprop.m] ================  


    
