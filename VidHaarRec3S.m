function Vid = VidHaarRec3S(DecVid,Lev)
% Directional Haar 3D Reconstruction: multi-Level
%
% Input
%    DecVid -  n1-by-n2-by-n3-by-m data
%    Lev    -  The decomposed level 
%
% Output
%    Vid    -  3d reconstructed data.
%                                                                       
% Copyright:
%    Xiaosheng Zhuang@CityU, 2018-Oct-9




if Lev == 1
    
   Vid              = HaarRec3S(DecVid,1);   
   
elseif Lev == 2     
    
   DecVid(:,:,:,6)  = HaarRec3S(DecVid(:,:,:,1:6), 2); 
   Vid              = HaarRec3S(DecVid(:,:,:,6:11),1); 
   
elseif Lev == 3    
    
   DecVid(:,:,:,6)  = HaarRec3S(DecVid(:,:,:,1:6), 3); 
   DecVid(:,:,:,11) = HaarRec3S(DecVid(:,:,:,6:11),2);   
   Vid              = HaarRec3S(DecVid(:,:,:,11:16),1);  
   
elseif Lev == 4
    
   DecVid(:,:,:,6) = HaarRec3S(DecVid(:,:,:,1:6), 4); 
   DecVid(:,:,:,11) = HaarRec3S(DecVid(:,:,:,6:11),3);   
   DecVid(:,:,:,16) = HaarRec3S(DecVid(:,:,:,11:16),2);   
   Vid              = HaarRec3S(DecVid(:,:,:,16:21),1);    
    
elseif Lev == 5

   DecVid(:,:,:,6)  = HaarRec3S(DecVid(:,:,:,1:6),  5);    
   DecVid(:,:,:,11) = HaarRec3S(DecVid(:,:,:,6:11), 4); 
   DecVid(:,:,:,16) = HaarRec3S(DecVid(:,:,:,11:16),3);   
   DecVid(:,:,:,21) = HaarRec3S(DecVid(:,:,:,16:21),2);   
   Vid              = HaarRec3S(DecVid(:,:,:,21:26),1);    

elseif Lev == 6
    
   DecVid(:,:,:,6)  = HaarRec3S(DecVid(:,:,:,1:6),  6);    
   DecVid(:,:,:,11) = HaarRec3S(DecVid(:,:,:,6:11), 5); 
   DecVid(:,:,:,16) = HaarRec3S(DecVid(:,:,:,11:16),4);   
   DecVid(:,:,:,21) = HaarRec3S(DecVid(:,:,:,16:21),3);   
   DecVid(:,:,:,26) = HaarRec3S(DecVid(:,:,:,21:26),2);
   Vid              = HaarRec3S(DecVid(:,:,:,26:31),1);     
    
end