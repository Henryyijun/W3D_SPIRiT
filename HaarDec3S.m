function DecVid = HaarDec3S(Vid,Lev)
% Directional Haar 3D Decomposition: 1 Level
%
% Input
%    Vid    -  Video 3D data. 
%    Lev    -  The decomposed level 
%
% Output
%    DecVid -  The result of video data.
%                                                                       
% Copyright:
%    Xiaosheng Zhuang@CityU, 2018-Oct-8
%

Vid0  = Vid;                                   % 0 0 0 

Vid1  = circshift(Vid0,[0,0,1]*2^(Lev-1));     % 0 0  1    -----> 1
Vid1m = circshift(Vid0,[0,0,-1]*2^(Lev-1));     % 0 0 -1    -----> -1

Vid2  = circshift(Vid0,[0,1,0]*2^(Lev-1));     % 0 1 0     -----> 2

Vid3  = circshift(Vid0,[0,1,1]*2^(Lev-1));     % 0  1  1   -----> 3
Vid3m = circshift(Vid0,[0,-1,-1]*2^(Lev-1));   % 0 -1 -1   -----> -3

Vid4  = circshift(Vid0,[1,0,0]*2^(Lev-1));     % 1 0 0     -----> 4

Vid5  = circshift(Vid0,[1,0,1]*2^(Lev-1));     %  1 0  1   ----->  5
Vid5m = circshift(Vid0,[-1,0,-1]*2^(Lev-1));   % -1 0 -1   -----> -5

Vid6  = circshift(Vid0,[1,1,0]*2^(Lev-1));         % 1 1 0     -----> 6

Vid7  = circshift(Vid0,[1,1,1]*2^(Lev-1));         % 1 1 1      -----> 7
Vid7m = circshift(Vid0,[-1,-1,-1]*2^(Lev-1));      % -1 -1 -1   -----> -7

Vidp1p1m1  = circshift(Vid0,[1,1,-1]*2^(Lev-1));   % 1 1 -1    -----> 
Vidm1m1p1  = circshift(Vid0,[-1,-1,1]*2^(Lev-1));  % -1 -1 1   -----> 

Vidp1m1p1  = circshift(Vid0,[1,-1,1]*2^(Lev-1));   % 1 -1 1     -----> 
Vidm1p1m1  = circshift(Vid0,[-1,1,-1]*2^(Lev-1));  % -1 1 -1   -----> 

Vidm1p1p1  = circshift(Vid0,[-1,1,1]*2^(Lev-1));   % -1 1 1     -----> 
Vidp1m1m1  = circshift(Vid0,[1,-1,-1]*2^(Lev-1));  % 1 -1 -1   -----> 

Vidp10m1  = circshift(Vid0,[1,0,-1]*2^(Lev-1));    % 1 0 -1   -----> 
Vidm10p1  = circshift(Vid0,[-1,0,1]*2^(Lev-1));    % -1 0 1   -----> 

Vid0p1m1  = circshift(Vid0,[0,1,-1]*2^(Lev-1));    % 0 1 -1   -----> 
Vid0m1p1  = circshift(Vid0,[0,-1,1]*2^(Lev-1));    % 0 -1 1   -----> 


% DecVid = zeros([size(Vid),6]);

DecVid(:,:,:,1)  = (Vid0+Vid1+Vid2+Vid3+Vid4+Vid5+Vid6+Vid7)/8; %low-pass filter

DecVid(:,:,:,2)  = (Vid4 - Vid0)/4; %(1,0,0)-(0,0,0) x-axis

DecVid(:,:,:,3)  = (Vid2 - Vid0)/4; %(0,1,0)-(0,0,0) y-axis

DecVid(:,:,:,4)  = sqrt(2)/8*(Vid6 - Vid0); %(1,1,0)-(0,0,0) xy-plane 1

DecVid(:,:,:,5)  = sqrt(2)/8*(Vid4 - Vid2); %(1,0,0)-(0,1,0) xy-plane 2

DecVid(:,:,:,6)  = 1/2*Vid0-1/16*(Vid1+Vid1m)-1/32*(Vid3+Vid3m)-...   % auxiliarly filter
                   1/32*(Vid5+Vid5m)-1/64*(Vid7+Vid7m)-1/64*(Vidp1p1m1+Vidm1m1p1)-...
                   1/64*(Vidp1m1p1+Vidm1p1m1)-1/64*(Vidm1p1p1+Vidp1m1m1)-1/32*(Vidp10m1+Vidm10p1)-...
                   1/32*(Vid0p1m1+Vid0m1p1);
return; 







