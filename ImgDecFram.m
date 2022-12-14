function [ResData]  =  ImgDecFram(Data, Frame)
%This function is to decompose the Image data by Frame operators.  
%
%Parameter:Data:  The image data to be decomposed. 
%                  Frame:  Frame operators. 
%                
%Return:      ResData:   The decomposed signal.             
%           
%Author:    galaxy 
%Date: 

DimFra  = size(Frame,3);     %How many operators 

for i = 1:DimFra
% ResData(:,:,i)  =  imfilter(Data,Frame(:,:,i),'replicate'); 
ResData(:,:,i)  =  imfilter(Data,Frame(:,:,i),'circular'); 
end

return; 