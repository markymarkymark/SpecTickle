function [hdr,path,file] = dicom_header(path,wildcard)
% Syntax: hdr = dicom_header(path,wildcard)
% Read Dicom headers of files matching wildcard
%
% To read all files in a folder matching a wildcard, try:
%	hdr=dicom_header('C:\temp\','*.dcm');
%
% If you only want to read a single file:
%	hdr=dicom_header('C:\temp\image.dcm');
%
% If you want to be prompted for a filename, make path = ''
%	e.g. hdr=dicom_header('');
%
% NOTE: the hdr structs are returned as a cell array, 
%    since MATLAB can't array different size structures
%
% Created: Mark A. Elliott, PhD
%   melliott@upenn.edu
%   https://www.med.upenn.edu/CAMIPM/mark-elliott.html

file = '';
if nargin < 1 , path     = ''	  ; end
if nargin < 2 , wildcard = '*.dcm'; end

% --- Prompt for user to choose a file, if none passed in ---
if isempty(path)
    filename = uigetfile_plus('dicom','dicompath',pwd(),{'*.dcm;*.MR;*.IMA';'*.*'},'Select a Dicom file');
    if (isempty(filename)), return; end
    hdr = dicominfo(filename, 'UseDictionaryVR',true);

% --- Read passed in filename ---
elseif (nargin < 2)
    if (~iscell(path))   % single file to read
        hdr = dicominfo(path, 'UseDictionaryVR',true); % This setting silences DicomDict warnings 
    else                    % cell array of filenames to read
        fullfiles = path;   % will read all files below
    end
end

% --- Done reading single dicom file? ---
if (exist('hdr'))
	return
end

% --- Returned will be a cell array ---
hdr    = {};

% --- User passed in cell array of complete filenames ---
if (exist('fullfiles')) 
    nfiles = numel(fullfiles);
    for i=1:nfiles
        hdrx = dicominfo(fullfiles{i}, 'UseDictionaryVR',true);
        hdr{i} = hdrx;	
    end
 
% --- Get filenames matching search args ---
else    
    files  = dir([path filesep() wildcard]);
    if isequal(wildcard,'*')	% remove '.' and '..' files returned if wildcard = '*'
        nfiles = numel(files);
        files = files(3:nfiles);
    end
    
    % --- Read the file(s) ---
    nfiles = numel(files);
    for i=1:nfiles
        hdrx = dicominfo([path filesep() files(i).name], 'UseDictionaryVR',true);
        hdr{i} = hdrx;	
    end
end
return