function [Ind]=roots_sorting(wave_number_imag,beta)

[row,col]=size(wave_number_imag);
nbeta=length(beta);

Ind=zeros(row,col);
for k=1:row
    Ind(k,1:end)=k;
end

% follow curves roots sorting
disp('roots sorting');
for k=3:nbeta-1
    for j=1:row
    p=polyfit(beta(k-2:k),[wave_number_imag(Ind(j,k-2),k-2),wave_number_imag(Ind(j,k-1),k-1),wave_number_imag(Ind(j,k),k)],2);
    F=polyval(p,beta(k+1));
    if (j>1)
        r=setdiff([1:row],Ind(1:j-1,k+1));
        [Y,I]=min(abs((F-wave_number_imag(r,k+1))));
        Ind(j,k+1)=r(I);
    else
        [Y,I]=min(abs((F-wave_number_imag(:,k+1))));
        Ind(j,k+1)=I;
    end
    
    end
end
