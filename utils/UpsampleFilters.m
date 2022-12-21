function [UpFilters]=UpsampleFilters(Filters,ZerosNum)
%This function is to upsample the filters by adding zeros. 
%For example: [1 2 3]->[1 0 2 0 3]
%
%Parameter:     Filters:   The filters band. 
%                     ZerosNum: The number of adding zeros.  
%Return:           UpFilters:  The result of upsampling the filters.
%                                                                       
%
%Auther:    galaxy 
%Date:       2008-10-28,(NUS).

% Filters = ones(3,3,2);
% Filters(3,3,2) = 2;
% ZerosNum = 1; 

[RowFil,ColFil,DimFil]= size(Filters);%The size of the filters. 
RowUpFilters = zeros(1,ColFil,DimFil);  
RowCount = 1; 
for i= 1: RowFil-1
  RowUpFilters(RowCount,:,:) = Filters(i,:,:); 
  RowUpFilters(RowCount+1:RowCount+ZerosNum,:,:) = 0; 
  RowCount = RowCount+ZerosNum+1;
end%End for 'i = 1: RowFil-1' 
RowUpFilters(RowCount,:,:) = Filters(i+1,:,:); 
%%%%
UpFilters = zeros(size(RowUpFilters,1),1,DimFil);  
ColCount = 1; 
for i= 1: ColFil-1
  UpFilters(:,ColCount,:) = RowUpFilters(:,i,:); 
  UpFilters(:,ColCount+1:ColCount+ZerosNum,:) = 0; 
  ColCount = ColCount+ZerosNum+1;
end%End for 'i = 1: RowFil-1' 
UpFilters(:,ColCount,:) = RowUpFilters(:,i+1,:); 
return; 






