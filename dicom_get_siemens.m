function [h,stat] = dicom_get_siemens(h)
% Get contents of Siemens private dicom MR header
% Private header fields are appended to passed in hdr
% returns:
%   0 if error
%   1 if Siemens private header succeccfully parsed
%   2 if no Siemens private header found
%
% Created: Mark A. Elliott, PhD
%   melliott@upenn.edu
%   https://www.med.upenn.edu/CAMIPM/mark-elliott.html

stat = 0;

% --- Check if Siemens scanner dicom ---
if (~isfield(h,'Manufacturer'))
    is_Siemens = 0;
else
    is_Siemens = ~isempty(findstr(lower(h.Manufacturer),'siemens'));
end
if (~is_Siemens)
    fprintf(1,'Dicom header is not from Siemens. Skipping private header parsing.\n');
    stat = 2;
    return
end

% --- Get names of all the Private fields ---
[status,val,h] = dicom_get_header(h);
if (status ~= 1)    % either an error (status=0) or no private header found status=2)
    stat = status;
    return 
end

% --- Run through all Siemens Private Image header fields ---
[np,slen] = size(h.PrivateImageNames);
for i=1:np
    tagname = deblank(h.PrivateImageNames(i,:));
    if (isfield(h,tagname))   % field already exists in public header
        %fprintf(1,'Ignoring duplicate fieldname "%s" in Public & Private header\n',tagname);
    else
        [status,val,h,dtype] = dicom_get_header(h,tagname);
        if (status == 0)
            return; 
        end
    end
end

% --- Run through all Siemens Private Sereis header fields ---
[np,slen] = size(h.PrivateSeriesNames);
for i=1:np
    tagname = deblank(h.PrivateSeriesNames(i,:));
    if (strncmp(tagname,'JUNK',4))
        disp('Ignoring unreadable PrivateSeries tag');
    elseif (isfield(h,tagname))   % field already exists in public header
        %fprintf(1,'Ignoring duplicate fieldname "%s" in Public & Private header\n',tagname);
    else
        [status,val,h,dtype] = dicom_get_header(h,tagname);
        if (status == 0)
            return; 
        end
    end
end

% --- extract the Sequence/Special card params (aka WipMemBlock[]) ---
if (isfield(h,'MrPhoenixProtocol'))
    max_doubles = 16;

    % --- wipmemblk DOUBLES ---
    Phoenix_lowercase = lower(h.MrPhoenixProtocol); % fix problem that siemens is not consistent on case!
%    h.WipMemBlock_Dvals = zeros(max_doubles,1);
    h.WipMemBlock_Dvals = cell(max_doubles,1);
    for i = 0:max_doubles-1
        %tagname     = sprintf('WipMemBlock_Dval_%02d',i);
        matchstring = sprintf('swipmemblock.adfree[%1d]',i);
        valstring   = [matchstring ' = %f'];
        val         = parse_phoenix(Phoenix_lowercase,matchstring, valstring ,1);
        if (isempty(val)), val = NaN'; end
        %h.(tagname)              = val;
        h.WipMemBlock_Dvals{i+1} = val; 
    end

    % --- wipmemblk LONGS ---
    max_longs   = 64;
%    h.WipMemBlock_Lvals = zeros(max_longs,1,'int32');
    h.WipMemBlock_Lvals = cell(max_longs,1);
    for i = 0:max_longs-1
        %tagname     = sprintf('WipMemBlock_Lval_%02d',i);
        matchstring = sprintf('swipmemblock.alfree[%1d]',i);
        valstring   = [matchstring ' = %d'];  % This MUST be %d to preserve int32 values
        val         = parse_phoenix(Phoenix_lowercase,matchstring, valstring ,1);
        if (isempty(val)), val = NaN; end
        %h.(tagname)              = val;
%        h.WipMemBlock_Lvals(i+1) = val; 
        h.WipMemBlock_Lvals{i+1} = val; 
    end

    % --- Did we encode param names in the unused LONG values? ---
    [h.WipMemBlock_Struct,h.WipMemBlock_LongLabels,h.WipMemBlock_DoubleLabels] = dicom_get_wipmemlabels(h.WipMemBlock_Lvals,h.WipMemBlock_Dvals);

    % --- get some other stuff ---
    h.('lSequenceID') = parse_phoenix(h.MrPhoenixProtocol,'lSequenceID',       'lSequenceID = %d'      ,1);
    h.('alTR0')       = parse_phoenix(h.MrPhoenixProtocol,'alTR[0]',           'alTR[0] = %d'          ,1);
    h.('alTR1')       = parse_phoenix(h.MrPhoenixProtocol,'alTR[1]',           'alTR[1] = %d'          ,1);  % for AFI sequence
    h.('SeqBinary')   = parse_phoenix(h.MrPhoenixProtocol,'tSequenceFileName', 'tSequenceFileName = %s',1);  
end

% --- Remove redundant fields created in "dicom_get_header() call ---
h = rmfield(h,'PrivateImageNames');
h = rmfield(h,'PrivateSeriesNames');
h = rmfield(h,'PrivateImageHdr');
h = rmfield(h,'PrivateSeriesHdr');

stat = 1;
end

% -------------------------------------------------------------------
function val = parse_phoenix(phoenixstring,matchstring,valstring,valindex)

val = [];
p = strfind(phoenixstring,matchstring);
if (~isempty(p))
    line = strtok(phoenixstring(p(1):p(1)+80),newline());
    c = textscan(line,valstring);
    val = c{valindex};
    if (iscell(val)), val = val{1}; end
    if (ischar(val))
        val = strrep(val,'"',''); 
        val = strrep(val,'%','%%'); 
    end
end
end



