%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%                       get_KM.m                      %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes global stiffness and mass matrix for a composite
% plate of given layup (layup) and material properties (C, rho)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [K, M] = get_KM_v2(nele_layer,np,h,C,rho)

nlayers = length(h);

% Number of elements per layer
nele = nele_layer * nlayers;

% 1D problem: mesh only in stacking (Z) direction
% DoF at each node = 3
nD = 1;                                     % do NOT change
ndof = 3;                                   % do NOT change

% Order of 1D elements
nnodes_ele = np+1;
nnodes = nele*(nnodes_ele-1) + 1;


%% IEN, LM and Coordinate arrays

% Nodal quadrature based on GLL points
[zi, w, SF, SFD] = GLLpoints(np);

% Populate IEN array (nD=1, hence dropping it)
IEN = zeros(nnodes_ele,nele);
LM = zeros(nnodes_ele*ndof,nele);
xcoord_IEN = zeros(nnodes_ele,nele); % actually z-coords

xcoord_IEN(1,1) = 0;
xcoord_IEN(end,end) = sum(h);

ele_no = 0;
node_no = 0;

for ii=1:nlayers
    layer_thick = h(ii);
    ele_length = layer_thick/nele_layer;
    
    for jj=1:nele_layer
        ele_no = ele_no + 1;
        xcoord_IEN(1,ele_no) = xcoord_IEN(1,1) + (ele_no-1)*ele_length;
        x_start = xcoord_IEN(1,ele_no);
        xcoord_IEN(nnodes_ele,ele_no) = xcoord_IEN(1,ele_no) + ele_length;
        x_end = xcoord_IEN(nnodes_ele,ele_no);
        
        for kk=1:(nnodes_ele-1)
            node_no = node_no + 1;
            if (kk>1)
                xcoord_IEN(kk,ele_no) = x_start * (1-zi(kk))/2 + x_end * (1+zi(kk))/2;
            end
            IEN(kk,ele_no) = node_no;
        end
        IEN(nnodes_ele,ele_no) = node_no+1;
    end
end

neq = nnodes * ndof;

for ii=1:nele
    for jj=1:nnodes_ele
        t1 = IEN(jj,ii);
        for kk=1:ndof
            t2 = (t1-1)*ndof + kk;
            t3 = (jj-1)*ndof + kk;
            LM(t3,ii) = t2;
        end
    end
end
      

%% Global K and M matrices

% Declare matrices
M = zeros(neq,neq);
K = cell(2,2);
for jj=1:2
    for kk=1:2
        K{jj,kk} = zeros(neq,neq);
    end
end

Ke = cell(2,2);

for i1=1:nlayers
    for i2=1:nele_layer
        
        % Declare element stiffness and mass matrices
        Me = zeros(size(LM,1),size(LM,1));
        for jj=1:2
            for kk=1:2
                Ke{jj,kk} = zeros(size(LM,1),size(LM,1));
            end
        end
        
        ele_no = (i1-1)*nele_layer + i2;
%         ele_length = xcoord_IEN(nnodes_ele,ele_no) - xcoord_IEN(1,ele_no);
        xcoord_ele = xcoord_IEN(:,ele_no);
                
        % Get Jacobian
        J = get_Jacobian(np,zi,SFD,xcoord_ele);
        
        % Get B and N matrices 
        for ii = 1:length(zi)
            [B, N] = get_B_v2(ii,np,J(ii),SF,SFD);
            for jj=1:2
                for kk=1:2
                    Ke{jj,kk} = Ke{jj,kk} + (B{jj})' * C{i1} * B{kk} * abs(J(ii)) * w(ii);
                end
            end
            
            Me = Me + rho * (N' * N) * abs(J(ii)) * w(ii);
        end
        
        % Assemble element stiffness and mass matrices
        M(LM(:,ele_no),LM(:,ele_no)) = M(LM(:,ele_no),LM(:,ele_no)) + Me;
        
        for jj=1:2
            for kk=1:2
                K{jj,kk}(LM(:,ele_no),LM(:,ele_no)) = K{jj,kk}(LM(:,ele_no),LM(:,ele_no)) + Ke{jj,kk};
            end
        end
        
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%