% Ong satin weave optimized constants
% W.H. Ong, N. Rajic, W.K. Chiu, C. Rosalie, Determination of the elastic properties of woven composite panels for
% Lamb wave studies, Composite Structures 141 (2016) 24–31
% Table 3

E11 = 44.17e9;
E22=E11;
E33 = 9.97e9;
G12 = 3.33e9;
G13 = G12;
G23 = G12;
ni12 = 0.35;
ni13 = ni12;
ni23 = 0.078;

ni21=ni12;
ni32=ni23;
G31=G13;

S= [1/E11     -ni21/E22        -ni21/E33     0         0         0 
     -ni12/E11   1/E22            -ni32/E33    0         0         0
     -ni13/E11   -ni23/E22       1/E33         0         0         0
       0            0                     0            1/G23    0        0
       0            0                     0             0       1/G31    0
       0            0                     0             0        0      1/G12]
 
C=inv(S)