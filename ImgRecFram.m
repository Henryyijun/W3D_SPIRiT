function [RecData]  =  ImgRecFram(Data, RecFrame)
%This function is to reconstruct the signal data by Frame operators.  
%
%Parameter: Data:  The sub-band image data to be reconstructed . 
%                  Frame:  Reconstuction Frame operators. 
%                
%Return:      RecData:   The reconstructed signal.             
%           
%Author:    galaxy 
%Date: 

DimFra  = size(RecFrame,3);     %How many operators 
RecData = zeros(size(Data(:,:,1)));

for i = 1:DimFra
%RecData  =  imfilter(Data(:,:,i),RecFrame(:,:,i),'symmetric') + RecData; 
RecData  =  imfilter(Data(:,:,i),RecFrame(:,:,i),'circular') + RecData; 
end

return; 