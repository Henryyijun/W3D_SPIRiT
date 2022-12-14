function DecVid = VidHaarDec3S(Vid,Lev)
% Directional Haar 3D Decomposition: multi Level
%
% Input
%    Vid    -  Video 3D data. 
%    Lev    -  The decomposed level 
%
% Output
%    DecVid -  The result of video data.
%                                                                       
% Copyright:
%    Xiaosheng Zhuang@CityU, 2018-Oct-9
%

if     Lev == 1  
    
    DecVid             = HaarDec3S(Vid,1); 
    
elseif Lev == 2
    
    DecVid2            = HaarDec3S(Vid,1);
    DecVid             = HaarDec3S(DecVid2(:,:,:,1),2);
    DecVid(:,:,:,7:11) = DecVid2(:,:,:,2:6);    
    
elseif Lev == 3
    
    DecVid3            = HaarDec3S(Vid,1);
    DecVid2            = HaarDec3S(DecVid3(:,:,:,1),2);
    DecVid             = HaarDec3S(DecVid2(:,:,:,1),3);
    DecVid(:,:,:,7:11) = DecVid2(:,:,:,2:6);
    DecVid(:,:,:,12:16)= DecVid3(:,:,:,2:6);
    
elseif Lev == 4
    
    DecVid4            = HaarDec3S(Vid,1);
    DecVid3            = HaarDec3S(DecVid4(:,:,:,1),2);
    DecVid2            = HaarDec3S(DecVid3(:,:,:,1),3);
    DecVid             = HaarDec3S(DecVid2(:,:,:,1),4);
    DecVid(:,:,:,7:11) = DecVid2(:,:,:,2:6);
    DecVid(:,:,:,12:16)= DecVid3(:,:,:,2:6); 
    DecVid(:,:,:,17:21)= DecVid4(:,:,:,2:6); 
    
elseif Lev == 5

    DecVid5            = HaarDec3S(Vid,1);    
    DecVid4            = HaarDec3S(DecVid5(:,:,:,1),2);
    DecVid3            = HaarDec3S(DecVid4(:,:,:,1),3);
    DecVid2            = HaarDec3S(DecVid3(:,:,:,1),4);
    DecVid             = HaarDec3S(DecVid2(:,:,:,1),5);
    DecVid(:,:,:,7:11) = DecVid2(:,:,:,2:6);
    DecVid(:,:,:,12:16)= DecVid3(:,:,:,2:6); 
    DecVid(:,:,:,17:21)= DecVid4(:,:,:,2:6);   
    DecVid(:,:,:,22:26)= DecVid5(:,:,:,2:6); 

elseif Lev == 6

    DecVid6            = HaarDec3S(Vid,1);
    DecVid5            = HaarDec3S(DecVid6(:,:,:,1),2);    
    DecVid4            = HaarDec3S(DecVid5(:,:,:,1),3);
    DecVid3            = HaarDec3S(DecVid4(:,:,:,1),4);
    DecVid2            = HaarDec3S(DecVid3(:,:,:,1),5);
    DecVid             = HaarDec3S(DecVid2(:,:,:,1),6);
    DecVid(:,:,:,7:11) = DecVid2(:,:,:,2:6);
    DecVid(:,:,:,12:16)= DecVid3(:,:,:,2:6); 
    DecVid(:,:,:,17:21)= DecVid4(:,:,:,2:6);   
    DecVid(:,:,:,22:26)= DecVid5(:,:,:,2:6); 
    DecVid(:,:,:,27:31)= DecVid6(:,:,:,2:6);
 
else
    error('Lev is too big');
end

return;





