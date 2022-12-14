function [ Cal_mat,Target_vec ] = Spirit_Kernel( Block,Acs )
%       Spirit_Kernel   creates the calibration matrix and the target vector
%                       for Spirit kernel
%       input:
%       Block           size of block
%       Acs             ACS data from multiple coils
%       output:
%       Cal_mat         calibration matrix
%       Target_vec      target vector
%
%       
%       Mach 18, 2018

%example
% Block = [3 3];
% % Acs = [1:6; 2:7; 3:8; 4:9; 5:10;  ];
% % Tem = 1:30;
% % Acs(:,:,2) =  reshape(Tem', [5,6]);
% Tem = 1:32;
% Acs =  reshape(Tem, [4,4,2]);
% clear Tem



[y,x,z] = size(Acs);

if(Block(1) > y || Block(2) > x)
    error('ACS area is too small!');
elseif (Block(1)<=1 || Block(2)<=1)
    error('block error!');
end

Cal_mat_row_length = (y-Block(1)+1)*(x-Block(2)+1); %How many blocks in each coil data 
Cal_mat_col_length = prod(Block);                               %The numbers of the parameters to be estimated   
Cal_mat = zeros(Cal_mat_row_length,Cal_mat_col_length,z);
count = 1;
for ii = 1:Block(2)                                                           %x direction  
        Block_x_head = ii;
        Block_x_end = x-Block(2)+ii;
    for jj = 1:Block(1)                                                      %y direction  
        Block_y_head = jj;
        Block_y_end = y-Block(1)+jj;
        Cal_mat(:,count,:) = reshape(Acs(Block_y_head:Block_y_end, Block_x_head:Block_x_end,:),[Cal_mat_row_length,z]);
        count = count + 1;
     end
end

Target_vec = reshape(Acs(ceil(Block(1)/2):y-floor(Block(1)/2), ceil(Block(2)/2):x-floor(Block(2)/2),:),[Cal_mat_row_length,z]);%Target point fiting 

%%
flag = 1;       
                
if(flag == 0)%0 means that the center piont is involved, to be computed 
  Cal_mat = reshape(Cal_mat,[Cal_mat_row_length,Cal_mat_col_length*z]);    
else%1 means that the center piont is not involved, not to be computed 
    Temp_mat =  reshape(Cal_mat,[Cal_mat_row_length,Cal_mat_col_length*z]); 
%     Cal_mat = cell(z,1);% reset Cal_mat to be a cell with size z(the number of coils)
    IndDel = floor(Block(1)-1)/2*Block(2)+ceil(Block(2)/2); %The index of the center point  of the first coil 
    Cal_mat = zeros(size(Temp_mat,1), size(Temp_mat,2)-1,z);
    
    for i = 1:z
        Target_loc = (i-1)*prod(Block)+IndDel;%delete location of the target point
        list = [1:Target_loc-1,Target_loc+1:prod(Block)*z];
        Cal_mat(:,:,i) = Temp_mat(:,list);
    end

end
%%
    

end

