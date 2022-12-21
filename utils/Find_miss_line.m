function [lines, new_missing_lines] = Find_miss_line(missing_lines, line_num)
% Function: Find_miss_line(missing_lines, line_num) find the lines that
%           need to be filled, extend from the middle
%
% Parameter: missing_lines: all the missing lines now, a coloumn
%            line_num: how many lines need to be finded
%
% Return: lines: which lines are finded
%         new_missing_lines: new missing lines update form missing_lines
%
% Sunrise
miss_size = size(missing_lines, 1);
if miss_size <= line_num
    lines = missing_lines;
    new_missing_lines = [];
else
    mid = floor(miss_size / 2);
    low = mid - floor(line_num / 2) + 1;
    high = mid + (line_num - floor(line_num / 2));
    lines = missing_lines(low: high);
    new_missing_lines = missing_lines;
    new_missing_lines(low: high) = [];
end
end
