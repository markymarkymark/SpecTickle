function [files,folderflag,finfo] = dir_plus(filesearch)
%  dir() func that returns full file or folder names, not structures
%  Also removes the '.' and '..' folders, if any
%
% Syntax:
%  [files,folderflag] = dir_plus(filesearch)
%
% Created: Mark A. Elliott, PhD
%   melliott@upenn.edu
%   https://www.med.upenn.edu/CAMIPM/mark-elliott.html
files      = {};
folderflag = [];
if (nargin < 1)
    finfo = dir();
else
    finfo = dir(filesearch);
end
if isempty(finfo), return; end
[files,folderflag,finfo] = finfo_plus(finfo);
end