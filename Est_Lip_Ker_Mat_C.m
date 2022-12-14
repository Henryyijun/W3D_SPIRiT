function  [Lip] = Est_Lip_Ker_Mat_C(Ful_kspace,Ker,Ker_Size)
%     Function:This  function is to Lipschitz estimation of the matrix  [ C] 
%       Input:
%       Ful_kspace: The full k-space data,the size is same as the FOV, with
%      the ACS data
%       Ker:   The interpolation kernel 
%       Ker_Size:    The kernel size as [7 ,7] 
%      
%       Output:
%      Lip: The lipschitz of the matrix C
%       June 9, 2018



Ker_Num = 1; 
 %%Generate the Fourier transfomer of kernel
[Row, Col, Coi] = size(Ful_kspace); 
Ker_Matrix = zeros(size(Ful_kspace,1),size(Ful_kspace,2)); 
for Coi_Num=1:Coi
for i=1:Coi
Ker_Matrix(1:Ker_Size(1), 1:Ker_Size(2)) = Ker(:,:,i,Coi_Num);
Ker_Matrix_Circ = circshift(Ker_Matrix, [-fix((Ker_Size(1)-1)/2) -fix((Ker_Size(2)-1)/2)]);
Ker_Matrix_Circ_FFT(:,:,Ker_Num) = fft2(Ker_Matrix_Circ)/sqrt(Row*Col);%;    
Ker_Num = Ker_Num+1; 
end
end


%%
% Ind = eye(Coi);
Eig_Val = zeros(Row*Col,1);
Eig_Num = 1;
for i=1:Row
    for j=1:Col
        Eig_Val(Eig_Num) =norm(reshape(Ker_Matrix_Circ_FFT(i,j,:),Coi,Coi));  
       
    end

end
Lip =  max(Eig_Val);

return;

















