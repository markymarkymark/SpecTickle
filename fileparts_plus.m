function [path,name,ext,fullname,botdir,topdir] = fileparts_plus(filename,force_cell)
% fileparts() function w/ added functionality
%
% Syntax:
%   [path,name,ext,fullname,botdir,topdir] = fileparts_plus(filename[,force_cell])
% Inputs:
%   fullname - complete filename (name + extension)
%   botdir   - bottom-most folder name
%   topdir   - the rest of the path above botdir
%
% Created: Mark A. Elliott, PhD
%   melliott@upenn.edu
%   https://www.med.upenn.edu/CAMIPM/mark-elliott.html

if (nargin < 2), force_cell = 0; end

[path,name,ext] = fileparts(filename);
fullname        = strcat(name,ext);
[topdir,botdir] = fileparts(path);

path   = strcat(path,filesep());
topdir = strcat(topdir,filesep());
botdir = strcat(botdir,filesep());

if (force_cell && ~iscell(path))
    path = {path};
    name = {name};
    ext  = {ext};
    fullname = {fullname};
    topdir   = {topdir};
    botdir   = {botdir};
end

end
