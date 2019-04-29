
Fabric_1
e22_f = 0.1*e11_f;
tic
[Q11,Q12,Q13,Q21,Q22,Q23,Q31,Q32,Q33,Q44,Q55,Q66,rho] = ...
    compfabricprop(fiberType, h_p, h_f, h_w, a_f, a_w, g_f, g_w, vol_0, ...
    e11_m, ni12_m, rho_m, e11_f, e22_f, ni12_f, ni23_f, rho_f,plot_weave);
toc