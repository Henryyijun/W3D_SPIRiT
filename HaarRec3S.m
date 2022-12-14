function Vid = HaarRec3S(DecVid,Lev)
% Directional Haar 3D Reconstruction: 1-Level
%
% Input
%    DecVid -  n1-by-n2-by-n3-by-14 data
%    Lev    -  The decomposed level 
%
% Output
%    Vid    -  3d reconstructed data.
%                                                                       
% Copyright:
%    Xiaosheng Zhuang@CityU, 2018-Oct-8

%% low-pass filter
Vid0  = DecVid(:,:,:,1);
Vid1  = circshift(Vid0,-[0,0,1]*2^(Lev-1));  % 0 0 1  -----> 1
Vid2  = circshift(Vid0,-[0,1,0]*2^(Lev-1));  % 0 1 0  -----> 2
Vid3  = circshift(Vid0,-[0,1,1]*2^(Lev-1));  % 0 1 1  -----> 3
Vid4  = circshift(Vid0,-[1,0,0]*2^(Lev-1));  % 1 0 0  -----> 4
Vid5  = circshift(Vid0,-[1,0,1]*2^(Lev-1));  % 1 0 1  -----> 5
Vid6  = circshift(Vid0,-[1,1,0]*2^(Lev-1));  % 1 1 0  -----> 6
Vid7  = circshift(Vid0,-[1,1,1]*2^(Lev-1));  % 1 1 1  -----> 7

Vid   = (Vid0+Vid1+Vid2+Vid3+Vid4+Vid5+Vid6+Vid7)/8;
%% (1,0,0)-(0,0,0) x-axis: 
Vid0  = DecVid(:,:,:,2);
Vid0  = circshift(Vid0,-[1,0,0]*2^(Lev-1))-Vid0;
Vid   = Vid + 1/4*Vid0;

%% (0,1,0)-(0,0,0) y-axis: 
Vid0  = DecVid(:,:,:,3);
Vid0  = circshift(Vid0,-[0,1,0]*2^(Lev-1))-Vid0;
Vid   = Vid + 1/4*Vid0;

%% (1,1,0)-(0,0,0) xy-plane 1:
Vid0  = DecVid(:,:,:,4);
Vid0  = circshift(Vid0,-[1,1,0]*2^(Lev-1))-Vid0;
Vid   = Vid + sqrt(2)/8*Vid0;

%% (1,0,0)-(0,1,0) xy-plane 2:
Vid0  = DecVid(:,:,:,5);
Vid0  = circshift(Vid0,-[1,0,0]*2^(Lev-1))-circshift(Vid0,-[0,1,0]*2^(Lev-1));
Vid   = Vid + sqrt(2)/8*Vid0;

%% auxiliarily filter 
Vid   = Vid + DecVid(:,:,:,6);
return; 







