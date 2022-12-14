function [NorCoe] = NormalizedCoeByEnergy(Coe);
%This is function to normalize the Fourier coefficients by energy.
%
%Parameter:      Coe:   The Fourier coefficients.                        
%                       
%Return:           NorCoe:   The normalized coefficients.     .
%                                                                
%
%Author:    galaxy 
%Date:      2013-3-11,(SZU)

EngCoe = abs(Coe).^(2);

SumEng = sum(sum(EngCoe,2),1);

%MaxEng = mean(SumEng(:))*10000000.25;
%MaxEng = mean(SumEng(:))*3.25;
MaxEng = (size(Coe,1)*size(Coe,2));

for i=1:size(SumEng,3)
   WeiEng(i) = MaxEng/SumEng(:,:,i);
end    
for i=1:size(SumEng,3)
   NorCoe(:,:,i) =Coe(:,:,i)*sqrt(WeiEng(i));
end 

% for i=1:size(SumEng,3)
%     Tem = Coe(:,:,i);
%     WeiEng(i) = abs(max(Tem(:)));
%     MinWei = abs(min(Tem(:)));
% end    
% for i=1:size(SumEng,3)
%    NorCoe(:,:,i) =Coe(:,:,i)/WeiEng(i);
% end 


 return;







