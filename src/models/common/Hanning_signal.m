function [t,st]=Hanning_signal(dt,nft,f_1,f_2,t_1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_t=1/f_1;   % total duration time of the excitation [s]
t_2=t_1+t_t; % excitation termiation time [s]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
st=zeros(nft,1);
t=zeros(nft,1);
for nf=1:nft
  st(nf)=0.0; t(nf)=(nf-1)*dt;
  if (t(nf) >= t_1) && (t(nf) <= t_2)
    st(nf)=0.5*(1-cos(2*pi*f_1*(t(nf)-t_2)))*sin(2*pi*f_2*(t(nf)-t_1));
    %st(nf)=0.5*(1-cos(2*pi*f_1*(t(nf)-t_2)))*cos(2*pi*f_2*(t(nf)-t_1)); %
    %cosine excitation
  end
end
