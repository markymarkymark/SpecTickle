function [path,foldernames,parentfolder,startpath] = uigetdir_plus(groupname,grouppath,defpath,varargin)
% uigetfile() function w/ added functionality like remembering the last folder selected
%   Also allows for selecting multiple folders by calling uigetdir2() (add 'multiselect','on')
%
% Syntax:
%   path = uigetdir_plus(groupname,grouppath,defpath,varargin...)
% Inputs:
%   groupname - the setpref() group to remember the last folder under (e.g. 'MyDicoms')
%   groupath  - the setpref() variable to remember (e.g. 'lastpath')
%   defpath   - the inital folder if none found under groupname/grouppath
%
% Created: Mark A. Elliott, PhD
%   melliott@upenn.edu
%   https://www.med.upenn.edu/CAMIPM/mark-elliott.html

path = ''; foldernames = ''; parentfolder = ''; startpath = '';
if (nargin < 2)
    fprintf(2,'USAGE: %s(groupname,grouppath,[defpath],varargs...)\n',mfilename())
    return
end
if (nargin < 3 || isempty(defpath)), defpath = pwd(); end
if (ispref(groupname,grouppath)), startpath = getpref(groupname,grouppath);
else,                             startpath = defpath; end

% Multiple folders requested? 
if (any(contains(lower(varargin),'multiselect'))) && (any(contains(lower(varargin),'on')))
    path = uigetdir2(startpath,varargin{1});
    if (isempty(path)), return; end % user hit cancel
else
    path = uigetdir(startpath,varargin{:});
    if (isequal(path,0)), path = ''; return; end % user hit cancel
end
%if (ischar(path)), path = {path}; end
[parentfolder,foldernames,ext] = fileparts_plus(path,1);
foldernames = strcat(foldernames,ext); % this handles '.' appearing in a folder name 
setpref(groupname,grouppath,parentfolder{1});   % save parent folder if multiple folders chosen
end