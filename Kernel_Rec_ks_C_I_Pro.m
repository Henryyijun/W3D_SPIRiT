function [ks_Rec]= Kernel_Rec_ks_C_I_Pro(ks_data, Ker,Pro)
%Function:This  function is to reconstuct the k-space data by the
%interplation kernel. 
%       Input:
%      ks_data: The full k-space data
%      Ker:   The interpolation kernel 
%      Pro:  The sampling positions of k-space       
%       Output:
%       ks_Rec: The  reconstructed each coil k-space data by  the   kernel.  
%       July 13, 2018

 for Coi = 1:size(ks_data,3)
      ks_Rec(:,:,Coi) = Pro.* (ImgRecFram(ks_data,Ker(:,:,:,Coi))-ks_data(:,:,Coi)); 
 end


return;