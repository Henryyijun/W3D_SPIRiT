function new_kspace = Fill_kspace(u, Sense, ks_data, Lines)
% Function: Fill_kspace(u, Sense, ks_data, Lines) fill some missing 
%           k-space lines (rows) of the center block
%
% Parameter: u: image
%            Sense: sensitivity, used to get k-space
%            ks_block: k-space block with missing lines, 3D
%            Lines: which k-space lines in ks_block are need to be filled
%            
% Return: new_kspace: new_kspace: filled k-space block
%
% Sunrise
current_ks_data = M_oper(u, Sense);
new_kspace = ks_data;
for i = 1: size(ks_data, 3)
    new_kspace(Lines, :, i) = current_ks_data(Lines, :, i);
end