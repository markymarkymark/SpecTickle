function [files,folderflag,finfo] = finfo_plus(finfo)
%  Takes finfo  struct from dir() and returns full file/flder names
%  Also removes the '.' and '..' folders, if any
%
% Syntax:
%  [files,folderflag] = finfo_plus(finfo)
%
% Created: Mark A. Elliott, PhD
%   melliott@upenn.edu
%   https://www.med.upenn.edu/CAMIPM/mark-elliott.html
files = [];
folderflag = [];
if (isempty(finfo)), return; end
nfiles     = numel(finfo);
pignore    = cellfun(@isequal,{finfo.name},repmat({'.'},1,nfiles)) | cellfun(@isequal,{finfo.name},repmat({'..'},1,nfiles)); % ignore '.' and '..'
files      = strcat({finfo(~pignore).folder},filesep(),{finfo(~pignore).name})';  % NOTE transpose so answer is nfiles x 1
folderflag = [finfo(~pignore).isdir]';
finfo      = finfo(~pignore);
end