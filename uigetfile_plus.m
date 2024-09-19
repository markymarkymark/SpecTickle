function [filename,path,name,ext,fullname,botdir,topdir] = uigetfile_plus(groupname,grouppath,defpath,varargin)
% uigetfile() function w/ added functionality like remembering the last folder selected
%
% Syntax:
%   [filename,path,name,ext,fullname,botdir,topdir] = uigetfile_plus(groupname,grouppath,defpath,varargin...)
% Inputs:
%   groupname - the setpref() group to remember the last folder under (e.g. 'MyDicoms')
%   groupath  - the setpref() variable to remember (e.g. 'lastpath')
%   defpath   - the inital folder if none found under groupname/grouppath
%
% Created: Mark A. Elliott, PhD
%   melliott@upenn.edu
%   https://www.med.upenn.edu/CAMIPM/mark-elliott.html

filename = ''; path = ''; name = ''; ext = ''; fullname = ''; botdir = ''; topdir = '';
if (nargin < 2)
    fprintf(2,'USAGE: %s(groupname,grouppath,[defpath],varargs...)\n',mfilename())
    return
end
if (nargin < 3), defpath = ''; end
if (isempty(defpath)), defpath = pwd(); end
filter = varargin{1};
prompt = varargin{2};

if (ispref(groupname,grouppath)), path = getpref(groupname,grouppath);
else, path = [defpath filesep()]; end

[file,path] = uigetfile(filter,prompt,path,varargin{3:end});
if (~iscell(file) && ~ischar(file)), path = ''; return; end   % user hit cancel
setpref(groupname,grouppath,path);                            % remember for next time

filename = strcat(path,file);
[path,name,ext,fullname,botdir,topdir] = fileparts_plus(filename,1);
end