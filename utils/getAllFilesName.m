function fileNames = getAllFilesName(path, str)
% 读取path 目录下所有以str结尾的文件
% 返回所有文件的文件名
% Author: Henry
% Date: 2021/4/25
% example 
% fileNames = getAllFilesName('F:\pMRI-3DTF', 'mat');
% fileNames =  {'brain_8ch.mat', 'phantom.mat', ...};
post = ['*.', str];
fileFolder = fullfile(path);
dirOutput = dir(fullfile(fileFolder, post));
fileNames = {dirOutput.name};
end