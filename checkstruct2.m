function [val,fieldname] = checkstruct2(s,varargin)
% Check a structure for one (or more) fieldname(s), with a default value to return if not found
% 
% Syntax:
%   value = checkstruct2(struct,fieldname1,[fieldname2],...,default_value)
%
% Created: Mark A. Elliott, PhD
%   melliott@upenn.edu
%   https://www.med.upenn.edu/CAMIPM/mark-elliott.html
%
val = [];
if (nargin < 3)
    fprintf(2,'ERROR:\n\tusage: "value = checkstruct2(struct,fieldname1,[fieldname2],...,default_value)"\n')
    return;
end
defval = varargin{end};
nfieldnames = nargin - 2;
for i=1:nfieldnames
    fieldname = varargin{i+1};
    if (isfield(s,fieldname))
        val = s.(fieldname);
        if (isempty(val)), val = defval; end  % empty vars get the default value
        return;
    end
end
fieldname = '';
val = defval;
end

