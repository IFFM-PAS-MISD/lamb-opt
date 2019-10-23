clear all; close all;

nlayers=8;
em_=zeros(50,1);
ef_=zeros(50,1);
Vol=zeros(50,1);
C11_=zeros(50,1);
ObjVal_=zeros(50,1);
for k_test=1:50
    
    load(['out',filesep,'test_',num2str(k_test),'_3.9mm_',num2str(nlayers),'lay_plain_wave.mat']);
    C11_(k_test) = C11;
    em_(k_test) =e11_m;
    ef_(k_test) =e11_f;
    Vol(k_test) =vol_0;
    ObjVal_(k_test)=ObjVal;
end

C11mean = mean(C11_)
[A,I]=min(ObjVal_)
%hist(em_,15);
%hist(ef_,15);
hist(Vol,15);
%hist(C11_,15);
%hist(ObjVal_,15)