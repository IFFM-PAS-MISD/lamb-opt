%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%                   get_Jacobian.m                    %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function gets Jacobian of element with given x coordinates
% (xcoord_ele) of given order (np). Integration point (also node) locations
% are given (zi) along with shape functions (SF) and their derivatives
% (SFD) at nodes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function J = get_Jacobian(np,zi,SFD,xcoord_ele)

intp = length(zi);
J = zeros(intp,1);

for inti=1:intp
    for nodei=1:np+1
        J(inti) = J(inti) + SFD(nodei,inti) * xcoord_ele(nodei);
    end
end 

end