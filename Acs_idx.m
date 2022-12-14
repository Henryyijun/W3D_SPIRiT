function [ Acs,AcsLoc ] = Acs_idx( KsData,NumOfAcsLine )
%   Acs_idx         create the ACS data and location
%   input:
%   KsData          k-space data
%   NumOfAcsLine    number of ACS lines
%   output:
%   Acs             ACS data
%   AcsLoc          location of ACS lines
%
%   Channing Lau
%   Dec. 2017

%Example
% KsData = rand(256,256,2);
% NumOfAcsLine = 32;
% NumOfAcsLine = [32 32];

[y,x,z] = size(KsData);
if length(NumOfAcsLine) == 1
    AcsLoc = floor((y-NumOfAcsLine+2)/2):floor((y+NumOfAcsLine)/2);
    Acs = KsData(AcsLoc,:,:);
elseif length(NumOfAcsLine) == 2
    AcsLoc(:,1) = floor((y-NumOfAcsLine(1)+2)/2):floor((y+NumOfAcsLine(1))/2);
    AcsLoc(:,2) = floor((x-NumOfAcsLine(2)+2)/2):floor((x+NumOfAcsLine(2))/2);
    Acs = KsData(AcsLoc(:,1),AcsLoc(:,2),:);
else
    error('illegal input');
end



end

