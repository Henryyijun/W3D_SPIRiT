function [Ker,Ker_Tra]=Kernel_Estimation(Ful_kspace,Ker_Size,ACS_Line)
%       Kernel_Estimation   creates the calibration matrix and the target vector
%                       for Spirit kernel
%       Input:
%      Ful_kspace: The full k-space data,the size is same as the FOV, with
%      the ACS data. 
%       Ker_Size:    The kernel size as [7 ,7] 
%       ACS_Line:  The number of lines around the center 
%       Output:
%       Ker:   The interpolation kernel 
%       Ker_Tra: The transpose kernel corresponding to Kernel matrix.  
%       May 5, 2018


% load ./data/05_t2_tse_tra_512_s33_3mm_8.mat
% Ful_kspace = ksfull;
% Ker_Size= [9 9];
%  ACS_Line = 32; 
 
 
 ACS_Ind = fix(size(Ful_kspace,1)/2) - fix(ACS_Line/2) + 1 :fix(size(Ful_kspace,1)/2) + fix(ACS_Line/2);
 ACS_Freq = 1:size(Ful_kspace,2);
 %ACS_Freq_Num = fix(size(Ful_kspace,2)/3);
 %ACS_Freq = fix(size(Ful_kspace,2)/2) - fix(ACS_Freq_Num/2) + 1 :fix(size(Ful_kspace,2)/2) + fix(ACS_Freq_Num/2);
  for i=1:size(Ful_kspace,3)
 ACS_Data(:,:,i)  =  Ful_kspace(ACS_Ind,ACS_Freq,i); 
 end
 %%
 [Mat_Data ,Tag_Vec ] = Spirit_Kernel( Ker_Size,ACS_Data);  % Construct the matrix data from the ACS data



 %%Estimate the coefficients in complex domain
%   for i=1:size(Mat_Data,3)
%  Mt = Mat_Data(:,:,i)';   
%  w(:,i) =  pinv(Mt*Mat_Data(:,:,i))*(Mt*Tag_Vec(:,i)); %Estimate the weight 
%   end
%%Estimate the coefficients in real domain 
 Mat_Data_r = real (Mat_Data);
 Mat_Data_i = imag (Mat_Data);
for i=1:size(Mat_Data,3)
 M =  [Mat_Data_r(:,:,i); Mat_Data_i(:,:,i)] ;
%  Mt = M';
 v  = [real(Tag_Vec(:,i)) ;imag(Tag_Vec(:,i)) ];
%  MtM = Mt*M;
%  Idn = eye(size(MtM,1));
%  w2(:,i) =  pinv(MtM+0.001*Idn)*(Mt*v);                                       %Estimate the weight 
  w(:,i)  = Solve_Linear_Constrain(M,v ,550);
%  Err(i) = norm(M*w(:,i) -v,2);
end
%  Err
 
Con_w = sum(w.^2,1);
Con_w_Max = max(Con_w);
% if(Con_w_Max>1)
% w = w/ sqrt(Con_w_Max);
% end
%%Construct the convolution kernel 
Ker_w = zeros(size(w,1)+1, size(w,2));
IndDel = floor(Ker_Size(1)-1)/2*Ker_Size(2)+ceil(Ker_Size(2)/2); %The index of the center point  of the first coil   
 for i = 1:size(Ker_w,2)
        Target_loc = (i-1)*prod(Ker_Size)+IndDel;                               %Deleted location of the target point
         Ind_List = [1:Target_loc-1,Target_loc+1:prod(Ker_Size)*size(Ker_w,2)];
        Ker_w(Ind_List,i) = w(:,i);
        
        Ker(:,:,:,i)=  reshape(Ker_w(:,i),[Ker_Size(1) Ker_Size(2) size(Ker_w,2)]);
        
        [RecFraOpe] =Transpose_Filter(Ker(:,:,:,i));
        Ker_Tra(:,:,:,i) =RecFraOpe ;
 end

 
 
 
  for coil = 1:size(Ker_w,2)
   for i= 1:size(Ker_w,2)     
       Trans_Ker(:,:,i,coil) =  Ker_Tra(:,:,coil,i) ;    
   end
  end
  
  Ker_Tra = Trans_Ker; %The transpose matrix of kernel matrix Ker;  
    
   KerW(:,1) = sum(sum(sum(Ker.^2,3),1),2) ;
   Ker_TraW(:,1)  = sum(sum(sum(Ker_Tra.^2,3),1),2);
   Ker_W_Mul = KerW.*Ker_TraW;
   Ker_W_Max =  max(Ker_W_Mul(:));
 if(Ker_W_Max>1)
 Ker_Con =sqrt(Ker_W_Max+eps);
 Ker = Ker/Ker_Con;
 Ker_Tra=Ker_Tra/Ker_Con;
 end

 
return;


function [xk] = Solve_Linear_Constrain(A,b,Itr)
% 
%
L= 1000;
xk = zeros(size(A,2),1);
y = xk;
AtA = A'*A;
Atb = A'*b;
tk =1; 

L = norm(AtA,Inf);

for i=1:Itr
    
   xk1 = y- (AtA*y-Atb)/L; 
   x_Square = norm(xk1,2); 
   if(x_Square>0.99)
       xk1= xk1/(x_Square*1.0102);     
   end
    %%%%%   
    tk1 = (1+sqrt(1+4*tk*tk))/2;
    Acc_Wei =(tk-1)/tk1;
    tk = tk1;   
   
   y = xk1+Acc_Wei*(xk1-xk); 
   xk = xk1; 
%    Err(i) = norm(A*xk -b,2);
end

return; 


function [RecFraOpe] =Transpose_Filter(FilOpe)
%This funciton is to transpose the filter 
% FilOpe: is 
RecFraOpe = FilOpe;
[RowFilOpe,ColFilOpe,DimFilOpe] = size(FilOpe);
HalfRow =  RowFilOpe/2;
HalfCol   =  ColFilOpe/2;
for i = 1:DimFilOpe    
         
       RowPlusPlus = round(HalfRow)+1;
       RowMinusMinus = fix(HalfRow) ;
       for j = 1:fix(HalfRow)  
             RecFraOpe(RowPlusPlus,:,i) = FilOpe(RowMinusMinus,:,i);   
             RecFraOpe(RowMinusMinus,:,i) = FilOpe(RowPlusPlus,:,i);   
             RowPlusPlus = RowPlusPlus+1;
             RowMinusMinus = RowMinusMinus -1 ;
       end   
end %Row sysmetric swap operation   
for i = 1:DimFilOpe   
       ColPlusPlus = round(HalfCol)+1;
       ColMinusMinus = fix(HalfCol) ;
       for j = 1:fix(HalfCol) 
             Temp = RecFraOpe(:,ColPlusPlus,i);
             RecFraOpe(:,ColPlusPlus,i) = RecFraOpe(:,ColMinusMinus,i);   
             RecFraOpe(:,ColMinusMinus,i) = Temp;   
             ColPlusPlus = ColPlusPlus+1;
             ColMinusMinus = ColMinusMinus -1 ;
       end   
end%Column sysmetric swap operation   

return;

















