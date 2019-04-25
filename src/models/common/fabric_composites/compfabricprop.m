function [Q11,Q12,Q13,Q21,Q22,Q23,Q31,Q32,Q33,Q44,Q55,Q66,rho] = ...
    compfabricprop(fiberType, h_p, h_f, h_w, a_f, a_w, g_f, g_w, vol_0, ...
    e11_m, ni12_m, rho_m, e11_f, e22_f, ni12_f, ni23_f, rho_f,plot_weave)
% COMPFABRICPROP   Calculation of the effective mechanical properties of fabric composite 
%    optional: more details about the function than in the H1 line 
%    optional: more details about the function than in the H1 line 
%    optional: more details about the function than in the H1 line 
% 
% Syntax: [output1,output2] = compfabricprop(input1,input2,input3) 
% 
% Inputs: 
%    fiberType - type of fabric, string: now available: 'plainWeave' 
%    h - ply thickness, dimensions [1, 1], Units: mm 
%    h_f - fill height, dimensions [1, 1], Units: mm
%    h_w - warp height, dimensions [1, 1], Units: mm
%    a_f - fill width, dimensions [1, 1], Units: mm
%    a_w - warp width, dimensions [1, 1], Units: mm
%    g_f - gap width along the fill, dimensions [1, 1], Units: mm
%    g_w - gap width along the warp, dimensions [1, 1], Units: mm
%    vol_0 - overall fiber volume fraction , dimensions [1, 1], Units: -
%    e11_m - matrix elastic modulus E11 , dimensions [1, 1], Units: GPa
%    ni12_m - matrix Poisson's ratio , dimensions [1, 1], Units: -
%    rho_m - matrix density , dimensions [1, 1], Units: kg/m3
%    e11_f - fiber elastic modulus E11 , dimensions [1, 1], Units: GPa
%    e22_f - fiber elastic modulus E22 , dimensions [1, 1], Units: GPa
%    ni12_f - fiber Poisson's ratio ni12, dimensions [1, 1], Units: -
%    ni23_f - fiber Poisson's ratio ni23, dimensions [1, 1], Units: -
%    rho_f - fiber density , dimensions [1, 1], Units: kg/m3
%    plot_weave - plot weave geometry if true, logical, troue or false, [-]
% 
% Outputs: 
%    Q11...Q66 - stiffness matrix components, integer, dimensions [m, n], Units: GPa 
%    rho - composite density, dimensions [1, 1], Units: kg/m3 
% Example: 
%    [output1,output2] = compfabricprop(input1,input2,input3) 
%    [output1,output2] = compfabricprop(input1,input2) 
%    [output1] = compfabricprop(input1,input2,input3) 
% 
% Other m-files required: plain_wave.m; plainweaveproperet.m; functionalHahn.m  
% Subfunctions: none 
% MAT-files required: none 
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2 
% 

% Author: Piotr Fiborek, D.Sc., Eng. 
% Institute of Fluid Flow Machinery Polish Academy of Sciences 
% Mechanics of Intelligent Structures Department 
% email address: pfiborek@imp.gda.pl 
% Website: https://www.imp.gda.pl/en/research-centres/o4/o4z1/people/ 
%%
%---------------------- BEGIN CODE----------------------

 h_p = h_p*1e-3;  h_f = h_f*1e-3;  h_w = h_w*1e-3; a_f = a_f*1e-3;  
 a_w = a_w*1e-3;  g_f = g_f*1e-3;  g_w = g_w*1e-3;
 rho = mean([vol_0*rho_f+(1-vol_0)*rho_m, (vol_0/rho_f+(1-vol_0)/rho_m)^-1]);
 %%
 Q = zeros(6);
 nNodes = 6;
 converCond = true;
 tol = 1e-2;
 while converCond
    switch fiberType
        case 'plainWeave'
         [H1,H2,H3,H4, Theta_f, Theta_w, vol_fw] = plain_weave(...
   h_p, h_f, h_w, a_f, a_w, g_f, g_w, vol_0,nNodes,plot_weave);  
        otherwise
            disp('')
    end
    % fiber volume fraction in [pure matrix, fill, warp, pure matrix]
    vol = [0 vol_fw vol_fw 0]; 
   
    [q11, q12, q13, q22, q23, q33, q44, q55, q66] = ...
    functionalHahn(e11_m, e11_f, e22_f, ni12_m, ni12_f, ni23_f, vol);
    
    
    % fill
    [Q11_f, Q12_f, Q13_f, Q21_f, Q22_f, Q23_f, Q31_f, Q32_f, Q33_f,...
        Q44_f, Q55_f, Q66_f] = ...
    plainweavepropert('fill', Theta_f, q11(:,2), q12(:,2), q13(:,2),...
    q22(:,2), q23(:,2), q33(:,2), q44(:,2), q55(:,2), q66(:,2));
    % warp
    [Q11_w, Q12_w, Q13_w, Q21_w, Q22_w, Q23_w, Q31_w, Q32_w, Q33_w,...
        Q44_w, Q55_w, Q66_w] = ...
    plainweavepropert('warp', Theta_w, q11(:,3), q12(:,3), q13(:,3),...
    q22(:,3), q23(:,3), q33(:,3), q44(:,3), q55(:,3), q66(:,3));
    % matrix
    [Q11_m, Q12_m, Q13_m, Q21_m, Q22_m, Q23_m, Q31_m, Q32_m, Q33_m,...
        Q44_m, Q55_m, Q66_m] = ...
    plainweavepropert('matrix', zeros(nNodes), q11(:,1), q12(:,1), q13(:,1),...
    q22(:,1), q23(:,1), q33(:,1), q44(:,1), q55(:,1), q66(:,1));
    %thicknesswise homogenization
    Q11_h = (Q11_m.*H1+Q11_f.*H2+Q11_w.*H3+Q11_m.*H4)/h_p;
    Q12_h = (Q12_m.*H1+Q12_f.*H2+Q12_w.*H3+Q12_m.*H4)/h_p;
    Q13_h = (Q13_m.*H1+Q13_f.*H2+Q13_w.*H3+Q13_m.*H4)/h_p;

    Q21_h = (Q21_m.*H1+Q21_f.*H2+Q21_w.*H3+Q21_m.*H4)/h_p;
    Q22_h = (Q22_m.*H1+Q22_f.*H2+Q22_w.*H3+Q22_m.*H4)/h_p;
    Q23_h = (Q23_m.*H1+Q23_f.*H2+Q23_w.*H3+Q23_m.*H4)/h_p;

    Q31_h = (Q31_m.*H1+Q31_f.*H2+Q31_w.*H3+Q31_m.*H4)/h_p;
    Q32_h = (Q32_m.*H1+Q32_f.*H2+Q32_w.*H3+Q32_m.*H4)/h_p;
    Q33_h = (Q33_m.*H1+Q33_f.*H2+Q33_w.*H3+Q33_m.*H4)/h_p;

    Q44_h = (Q44_m.*H1+Q44_f.*H2+Q44_w.*H3+Q44_m.*H4)/h_p;
    Q55_h = (Q55_m.*H1+Q55_f.*H2+Q55_w.*H3+Q55_m.*H4)/h_p;
    Q66_h = (Q66_m.*H1+Q66_f.*H2+Q66_w.*H3+Q66_m.*H4)/h_p;
    
   Q11 = sum(sum(Q11_h))/numel(Q11_h);
   Q12 = sum(sum(Q12_h))/numel(Q12_h);
   Q13 = sum(sum(Q13_h))/numel(Q13_h);
   
   Q21 = sum(sum(Q21_h))/numel(Q21_h);
   Q22 = sum(sum(Q22_h))/numel(Q22_h);
   Q23 = sum(sum(Q23_h))/numel(Q23_h);
   
   Q31 = sum(sum(Q31_h))/numel(Q31_h);
   Q32 = sum(sum(Q32_h))/numel(Q32_h);
   Q33 = sum(sum(Q33_h))/numel(Q33_h);
   
   Q44 = sum(sum(Q44_h))/numel(Q44_h);
   Q55 = sum(sum(Q55_h))/numel(Q55_h);
   Q66 = sum(sum(Q66_h))/numel(Q66_h);
  % convergence condition 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  iQ = [Q11 Q12 Q13 0 0 0
      Q21 Q22 Q23 0 0 0
      Q31 Q23 Q33 0 0 0
      0 0 0 Q44 0 0
      0 0 0 0 Q55 0
      0 0 0 0 0 Q66];
  
    if any(any(abs(Q-iQ)>tol))
       Q = iQ;
       nNodes = nNodes+4;
    else
        converCond = false;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
 end

%---------------------- END OF CODE---------------------- 

% ================ [compfabricprop.m] ================  
