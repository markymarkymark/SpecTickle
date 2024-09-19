function [stat,val,hdr,dtype] = dicom_get_header(hdr,fieldname)
% Return a Dicom header value
%	Checks for field in standard header and Siemens Private header.
%   If called w/o a fieldname, then extracts ALL Siemens Private header fields
% returns:
%   0 if error
%   1 if Siemens private header succeccfully parsed
%   2 if no Siemens private header found
%
% Created: Mark A. Elliott, PhD
%   melliott@upenn.edu
%   https://www.med.upenn.edu/CAMIPM/mark-elliott.html

stat  = 0;
val   = 0;
dtype = '';

% -------------------------------------------------
% --- Return existing header field if it exists ---
% -------------------------------------------------
if (nargin > 1)
    if (isfield(hdr,fieldname))
        val = getfield(hdr,fieldname);
        stat = 1;
        return;
    end
end

% -------------------------------------------------
% --- Look for field in Siemens Private header ---
% -------------------------------------------------
%	fprintf(1,'No match for "%s" in standard Dicom header fields. Searching Private headers...\n',fieldname);

% --- Extract the Siemens Private headers ---
if (~isfield(hdr,'Private_headers_parsed'))
		
    % --- is this an MRS file? ---
    hdr.is_spec = isfield(hdr,'Private_7fe1_1010');
    
    % --- is this a CSI file ? ---
    if (~isfield(hdr,'Rows')), hdr.Rows = 1; end  % Need to fix this for CSI
    hdr.is_csi = (hdr.is_spec & (hdr.Rows > 1));

    % --- Call SPM module to parse Private fields ---
    tmp  = me_spm8_dicom_headers(hdr.Filename);
    if (isempty(tmp)), return; end
    thdr = tmp{1};
		
    % --- find and rename the Private Image header fields ---
    if (isfield(thdr,'Private_0029_1010'))
        hdr.PrivateImageHdr = thdr.Private_0029_1010;
    elseif (isfield(thdr,'Private_0029_1110'))
        hdr.PrivateImageHdr = thdr.Private_0029_1110;
    elseif (isfield(thdr,'Private_0029_1210'))
        hdr.PrivateImageHdr = thdr.Private_0029_1210;
    elseif (isfield(thdr,'CSAImageHeaderInfo'))
        hdr.PrivateImageHdr = thdr.CSAImageHeaderInfo;
    elseif (isfield(thdr,'CSANonImageHeaderInfoVB'))
        hdr.PrivateImageHdr = thdr.CSANonImageHeaderInfoVB;
    end

    % --- find and rename the Series Image header fields ---
    if     (isfield(thdr,'Private_0029_1020'))   && (numel(thdr.Private_0029_1020) > 1)    hdr.PrivateSeriesHdr = thdr.Private_0029_1020;
    elseif (isfield(thdr,'Private_0029_1120'))   && (numel(thdr.Private_0029_1120) > 1)    hdr.PrivateSeriesHdr = thdr.Private_0029_1120;
    elseif (isfield(thdr,'Private_0029_1220'))   && (numel(thdr.Private_0029_1220) > 1)    hdr.PrivateSeriesHdr = thdr.Private_0029_1220;
    elseif (isfield(thdr,'CSASeriesHeaderInfo')) && (numel(thdr.CSASeriesHeaderInfo) > 1)  hdr.PrivateSeriesHdr = thdr.CSASeriesHeaderInfo;
    elseif (isfield(thdr,'CSAMiscProtocolHeaderInfoVB')) && (numel(thdr.CSAMiscProtocolHeaderInfoVB) > 1)  hdr.PrivateSeriesHdr = thdr.CSAMiscProtocolHeaderInfoVB;
    end

    % --- Error if didn't find Private headers ---
    if (~isfield(hdr,'PrivateImageHdr') || ~isfield(hdr,'PrivateSeriesHdr'))
        fprintf(1,'Did not find find private Image or private Series headers.\n');
        stat = 2;
        return;
    end

    % --- Remember so we don't need to do this again ---
    hdr.Private_headers_parsed = 1;
    if (isfield(hdr,'PrivateImageHdr')),  hdr.PrivateImageNames  = strvcat(hdr.PrivateImageHdr.name); end
    if (isfield(hdr,'PrivateSeriesHdr')), hdr.PrivateSeriesNames = strvcat(hdr.PrivateSeriesHdr.name); end
end

% --- Return if no fieldname asked for ---
if (nargin < 2)
    stat = 1;
    return
end

% --- Search Private headers for desired field ---
index = [];
if (isfield(hdr,'PrivateImageHdr'))
    index = strmatch(fieldname,strvcat(hdr.PrivateImageHdr.name),'exact');
    if (~isempty(index)), entry = hdr.PrivateImageHdr(index); end
end
if ((isfield(hdr,'PrivateSeriesHdr')) && isempty(index))
    index = strmatch(fieldname,strvcat(hdr.PrivateSeriesHdr.name),'exact');
    if (~isempty(index)), entry = hdr.PrivateSeriesHdr(index); end
end
	
% --- No match ---
if (isempty(index))
    fprintf(2,'Did not find match for "%s" in Private Dicom header fields.\n',fieldname);
    return;

% --- empty field value ?? ---
elseif (~isfield(entry,'nitems'))
    %%fprintf(1,'Found match for "%s" in Private Dicom header fields cant read this format.\n',fieldname);
    return;
    
elseif (entry.nitems == 0)
    %%fprintf(1,'Found match for "%s" in Private Dicom header fields but value is missing.\n',fieldname);
    val = '';
    dtype = '%s';

% --- Parse field value(s) ---
else
    if (entry.vm == 0), entry.vm = entry.nitems; end		% BUG!!!

    % --- Pull out value(s) ---
    switch(deblank(entry.vr))
        case {'DS','IS','US','FD','UI','UL','SL'},
            for i=1:entry.vm, val(i) = str2double(entry.item(i).val); end
            val = val';
            dtype = '%g';
        case {'AE','AS','AT','CS','DA','DS','DT','FL','LO','LT','PN','SH','SS','ST','TM','UT'},
            val = deblank(strvcat(entry.item(1:entry.vm).val));
            dtype = '%s';
        case {'UN'},    % This is the SIEMENS "ASCCONV" field!!!!!
            val = strvcat(entry.item(1:entry.vm).val);
            dtype = '%20s';
        otherwise,
            fprintf(2,'Unknown fieldtype %s.\n',entry.vr);
            return;
    end
end

% --- Save result and done ---
fieldname = strrep(fieldname,'-', '_'); % fix problem if fieldname has a dash in it.
hdr = setfield(hdr,fieldname,val);
stat = 1;

return;

