%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%                     GLLpoints.m                     %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes GLL node locations (zi), corresponding shape
% functions (SF), their derivatives (SFD), weight of nodal quadrature (w)
% for a given order of sprectral elements (np)
% SF(ii,jj) = value of ii_th shape function at zi_jj
% SFD(ii,jj) = derivative of ii_th shape function at zi_jj
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [zi, w, SF, SFD] = GLLpoints(np)


%% 3 <= np <= 5

if (np < 3)
    np = 3;
    disp('taking np = 3');
elseif (np > 5)
    np = 5;
    disp('taking np = 5');
end


%% Declare matrices

zi = zeros(np+1,1);
w = zeros(np+1,1);
SF = ones(np+1,np+1);
SFD = zeros(np+1,np+1);


%% GLL quadrature table

if (np==3)
    zi(1) = -1.0;
    zi(2) = -1/5*sqrt(5);
    zi(3) = 1/5*sqrt(5);
    zi(4) = 1;
    
    w(1) = 1/6;
    w(2) = 5/6;
    w(3) = 5/6;
    w(4) = 1/6;
    
elseif (np==4)
    zi(1) = -1.0;
    zi(2) = -1/7*sqrt(21);
    zi(3) = 0;
    zi(4) = 1/7*sqrt(21);
    zi(5) = 1.0;
    
    w(1) = 1/10;
    w(2) = 49/90;
    w(3) = 32/45;
    w(4) = 49/90;
    w(5) = 1/10;
    
elseif (np==5)
    zi(1) = -1.0;
    zi(2) = -sqrt(1/21*(7+2*sqrt(7)));
    zi(3) = -sqrt(1/21*(7-2*sqrt(7)));
    zi(4) = sqrt(1/21*(7-2*sqrt(7)));
    zi(5) = sqrt(1/21*(7+2*sqrt(7)));
    zi(6) = 1.0;
    
    w(1) = 1/15;
    w(2) = 1/30*(14-sqrt(7));
    w(3) = 1/30*(14+sqrt(7));
    w(4) = 1/30*(14+sqrt(7));
    w(5) = 1/30*(14-sqrt(7));
    w(6) = 1/15;
end


%% Shape functions and their derivatives

denom = ones(np+1,1);

for ii=1:np+1
    for jj=1:np+1
        if(jj~=ii)
            denom(ii) = denom(ii) * (zi(ii) - zi(jj));
        end
    end
end

for ip=1:length(zi)
    for ii=1:np+1
        for jj=1:np+1
            if(jj~=ii)
                SF(ii,ip) = SF(ii,ip) * (zi(ip) - zi(jj));
            end
        end
        SF(ii,ip) = SF(ii,ip) / denom(ii);
    end
end

for ip=1:length(zi)
    for ii=1:np+1
        t1 = 0;
        for jj=1:np+1
            if(jj~=ii)
                t2 = 1;
                for kk=1:np+1
                    if((kk~=ii) && (kk~=jj))
                        t2 = t2 * (zi(ip) - zi(kk));
                    end
                end
                t1 = t1 + t2;
            end
        end
        SFD(ii,ip) = t1/denom(ii);
    end
end


%% Verifying with analytical
% 
% syms x;
% 
% f = cell(np+1,1);
% 
% for ii=1:np+1
%     f{ii} = 1;
% end
% 
% for ii=1:np+1
%     for jj=1:np+1
%         if(jj~=ii)
%             f{ii} = f{ii}*(x-zi(jj));
%         end
%     end
%     f{ii} = f{ii} / denom(ii);
% end
% 
% for ii=1:np+1
%     for jj=1:np+1
%         SFA(ii,jj) = subs(f{ii},zi(jj));
%         SFDA(ii,jj) = subs(diff(f{ii},x),zi(jj));
%     end
% end
% 
% SFA
% SFDA       


end   

    
    