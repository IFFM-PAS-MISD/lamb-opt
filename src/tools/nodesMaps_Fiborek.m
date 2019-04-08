function [I_G,I_L]=nodesMaps_Fiborek(elementNodes,n,n_z)
numberNodes = max(max(elementNodes));
numberElements = size(elementNodes,1);
Psi = zeros(numberNodes,1);
phi_E = zeros(numberNodes,12);
phi_L = zeros(numberNodes,12);

for e = 1:numberElements
    for i = 1:n^2*n_z;
        k = elementNodes(e,i);
        Psi(k) = Psi(k)+1;
        phi_E(k,Psi(k)) = e;
        phi_L(k,Psi(k)) = i;
    end
end
I_G = zeros(n^2*n_z*numberElements/12,12);
I_L = zeros(n^2*n_z*numberElements/12,12);
H = zeros(12,1);
for k = 1:numberNodes
   [~,S] = sort(H);
   for i = 1:Psi(k)
       H(S(i)) = H(S(i))+1;
       e = phi_E(k,i);
       I_G(H(S(i)),S(i)) = k;
       I_L(H(S(i)),S(i)) = (e-1)*n^2*n_z+phi_L(k,i);
   end
end
  
    
