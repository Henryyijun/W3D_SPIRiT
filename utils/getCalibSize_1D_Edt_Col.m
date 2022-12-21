function [calibSize, ACS_edge]= getCalibSize_1D_Edt_Col(Mask)
%   calibSize = getCalibSize_edt(mask)
%
%   Givem a 1D sampling mask, the function attempts to estimate the largest region in 
%   the center of k-space that is fully sampled in order to use for calibraion.
%   Simple version of getCalibSize by Michael Lustig 2009.
%
%   The Mask is sampled in column, different from the function
%   getCalibSize_1D_Edt.m
%   
% Input:
%		mask - 2D binary array
%
% Output:	
%		calibSize - [1x2] size of the calibration area
%

mask_row = size(Mask,1);
mask_col = size(Mask,2);
Mask_row_center = mask_row/2;
Mask_col_center = mask_col/2;
sx = 1;
sy = 1;
% row
row_pos = Mask_row_center - 1;
while row_pos > 0
    if Mask(row_pos, Mask_col_center) == 0
        break
    else
        sx = sx + 1;
    end
    row_pos = row_pos - 1;
end
ACS_row_up_edge = row_pos + 1;
row_pos = Mask_row_center + 1;
while row_pos <= mask_row
    if Mask(row_pos, Mask_col_center) == 0
        break
    else
        sx = sx + 1;
    end
    row_pos = row_pos + 1;
end
ACS_row_down_edge = row_pos - 1;
% col
col_pos = Mask_col_center - 1;
while col_pos > 0
    if Mask(Mask_row_center, col_pos) == 0
        break
    else
        sy = sy + 1;
    end
    col_pos = col_pos - 1;
end
ACS_col_left_edge = col_pos + 1;
col_pos = Mask_col_center + 1;
while col_pos <= mask_col
    if Mask(Mask_row_center, col_pos) == 0
        break
    else
        sy = sy + 1;
    end
    col_pos = col_pos + 1;
end
ACS_col_right_edge = col_pos - 1;
calibSize = [sx,sy];
ACS_edge = [ACS_col_left_edge, ACS_col_right_edge];