function par = svs_acqparams(files,opt)
% Get structure of important SVS params, including DFS specific
% "files" can be:
%   Dicom filename
%   Dicom header (from dicomread())
%   TWIX data filename
%   TWIX hdr object (from mapvbvd())
%       if filename(s), it should be a cell array of filenames

par = [];
if (nargin < 1), files     = {}; end
if (nargin < 2), opt.dummy = 1;  end
opt.verbose         = checkstruct(opt,'verbose',1);
opt.m0_duration     = checkstruct(opt,'m0_duration',5000);     % Assign this TI/TS time for M0 scans
opt.remove_oversamp = checkstruct(opt,'remove_oversamp',0);    % remove OS

if (isempty(files))
    files = uigetfile_plus(mfilename(),'lastpath',pwd(),{'*.dcm;*.MR;*.IMA;*.dat';'*.*'},'Select a Dicom or TWIX file');
    if (isempty(files)), return; end
end
if (~iscell(files)), files = {files}; end % make input a cell array of filenames/structs/objects

% --- Dicom or TWIX raw data? ---
exampfile = files{1};
switch (class(exampfile))
    case 'char'
        [~,~,ext] = fileparts(exampfile);
        switch lower(ext)
            case '.dat'
                is_dcm = 0;
            case {'.dcm','.mr','.ima'}
                is_dcm = 1;
            otherwise
                fprintf(2,'ERROR: unrecognized filename extension type %s\n',ext);
                return
        end
    case 'struct'
        if (isfield(exampfile,'hdr'))
            is_dcm = 0;
        else
            is_dcm = 1;
        end
    otherwise
        fprintf(2,'ERROR: unrecognized input type %s\n',class(exampfile));
        return
end

% --- Get params ---
if (is_dcm), par = dcm_acqparams(files, opt);
else,        par = raw_acqparams(files, opt); end

end

% ------------------------------------------------------------------------------------------------
% Get params from DICOM file
% ------------------------------------------------------------------------------------------------
function par = dcm_acqparams(dcmfiles, opt)

% -- make structure template ---
nfiles = numel(dcmfiles);
par    = make_par_struct(nfiles);

% --- now fill struct array ---
for i=1:nfiles
    if (ischar(dcmfiles{i}))
        h = dicom_header(dcmfiles{i});
    else
        h = dcmfiles{i};
    end
    if (~isfield(h,'Private_headers_parsed')), h = dicom_get_siemens(h); end  % this reads private headers AND MrPhoenixProtocol

    par(i).is_EJA           = contains(lower(h.SequenceName),'eja');
    par(i).is_DFS           = contains(lower(h.SeqBinary),'svssel') || contains(lower(h.SeqBinary),'s3sel') || contains(lower(h.SeqBinary),'svs_metab') || contains(lower(h.SeqBinary),'selspec');

    par(i).TR               = h.RepetitionTime;
    par(i).TE               = h.EchoTime;
    par(i).nt               = h.DataPointColumns;
    par(i).nex              = h.NumberOfAverages;
    par(i).actual_nex       = h.NumberOfAverages;
    par(i).dwelltime        = h.RealDwellTime;
    par(i).BW               = 1/(par(i).dwelltime * 1e-9);  % NOTE that we assume here that removeOS was OFF
    par(i).f0               = h.ImagingFrequency;
    par(i).txrefvolt        = checkstruct(h,'TransmitterReferenceAmplitude',h.TransmitterCalibration);
    par(i).seriesnum        = h.SeriesNumber;
    par(i).imagenum         = h.AcquisitionNumber;
    par(i).seriesname       = string(h.SeriesDescription);
    par(i).voxdims          = [h.VoiThickness, h.VoiPhaseFoV, h.VoiReadoutFoV];
    par(i).voxsize          = prod(par(i).voxdims)/1000;                     % voxel size (cc)
    par(i).scandate         = str2num(h.StudyDate);
    %    par(i).scantime        = fix(str2num(h.SeriesTime));
    par(i).scantime         = fix(str2num(h.StudyTime));  % same for all series
    par(i).patname          = h.PatientName.FamilyName; if (isempty(par(i).patname)), par(i).patname = 'xxx'; end
    par(i).coilname         = h.TransmittingCoil;
    par(i).nucleus          = h.ImagedNucleus;
    par(i).fieldstrength    = [num2str(h.MagneticFieldStrength) 'T'];
    par(i).seqbinary        = h.SeqBinary;
    par(i).scanID           = h.lSequenceID; % this is tied to the TWIX raw data "MID"
    if (isfield(h,'PatientAge'))
        par(i).patage = sscanf(h.PatientAge,'%dY');
    else
        par(i).patage = str2num(h.StudyDate(1:4)) - str2num(h.PatientBirthDate(1:4)); % unbelievable!!
    end

    % --- Figure out DFS sequence info from Special card params ---
    if (par(i).is_DFS)

        par(i) = get_dfs_params(par(i), h.WipMemBlock_Struct, h.WipMemBlock_Lvals, h.WipMemBlock_Dvals, opt);

        % --- Did we encode wipmem params with labels? ---
        % if (~isempty(h.WipMemBlock_Struct))
        %     seqtypes = {'Voxel','Slice'};
        %     rfpulses = {'eburp','hsinc',sprintf('UserSinc%1dx%1d',h.WipMemBlock_Struct.SncL,h.WipMemBlock_Struct.SncR)};
        %     explist  = {'svs', 'inv_selsinc', 'inv_selHS', 'inv_nonsel', 'sat_sel', 'sat_nonselAHP', 'sat_nonselBIR4', 'sat_nonselRECT', 'sat_nonselSINC'};
        %     par(i).seqver       = h.WipMemBlock_Struct.Vers;
        %     par(i).presamp      = h.WipMemBlock_Struct.DumP;
        %     par(i).excitepulse  = rfpulses{h.WipMemBlock_Struct.ExPu+1};
        %     par(i).seqtype      = seqtypes{h.WipMemBlock_Struct.LocM+1};
        %     par(i).exptype      = explist{h.WipMemBlock_Struct.InvM+1};
        %     par(i).presamp      = h.WipMemBlock_Struct.DumP * 2;
        %     par(i).adcdelay     = h.WipMemBlock_Struct.ADCD;
        %     par(i).exciteppm    = h.WipMemBlock_Struct.ExCn;
        %     par(i).excitewidth  = h.WipMemBlock_Struct.ExWd;
        %     par(i).TI           = h.WipMemBlock_Struct.TSTI;  if (h.WipMemBlock_Struct.InvM ==  0), par(i).TI = opt.m0_duration; end
        %     par(i).crushcycle   = h.WipMemBlock_Struct.CyCr;
        %
        % % --- Have to guess what the params are ---
        % else
        %     h.WipMemBlock_Lvals = double(h.WipMemBlock_Lvals);
        %
        %     %version = h.WipMemBlock_Lvals(1);
        %     version = h.WipMemBlock_Dvals(1);       % newer versions wrote the version number to Dval00
        %     if (version > 0 && version < 1)
        %         par(i).seqver = version + 10;
        %     else
        %         lval00 = h.WipMemBlock_Lvals(1);     % Klugey way to figure out SVSSEL version from older versions
        %         switch(lval00)
        %             case 2
        %                 par(i).seqver = 1.0;
        %             case -1
        %                 par(i).seqver = 2.0;
        %             case 1
        %                 dval07 = h.WipMemBlock_Dvals(8);
        %                 lval05 = h.WipMemBlock_Lvals(6);
        %                 if (dval07 ~= -1),     par(i).seqver = 3.0;         % OLD voxel-selective w/ samples-before-echo
        %                 elseif (lval05 == -1), par(i).seqver = 2.0;
        %                 else,                  par(i).seqver = 5.0; end     % OLD slice-selective w/ samples-before-echo
        %             otherwise
        %                 fprintf(2,'ERROR: Cant figure sequence version from h.WipMemBlock_Lvals(1) = %1d\n',lval00);
        %                 return
        %         end
        %     end
        %
        %     % --- Parse based on sequence version ---
        %     switch(par(i).seqver)
        %         case 1.0
        %             par(i).excitepulse = 'eburp';
        %             par(i).refocpulse  = 'lobw';
        %             par(i).seqtype     = 'Voxel';
        %             par(i).TI          = h.WipMemBlock_Lvals(5);  % This is ALSO the Sat. Recovery time!!
        %             par(i).crusher     = h.WipMemBlock_Dvals(3);
        %             par(i).exciteppm   = h.WipMemBlock_Dvals(4);  % center of excite pulse (ppm)
        %             par(i).excitewidth = h.WipMemBlock_Dvals(5);  % width of excite pulse (ppm)
        %         case 2.0
        %             if (lval00 == -1), par(i).excitepulse = 'eburp';
        %             else               par(i).excitepulse = 'hsinc'; end
        %             lval05 = h.WipMemBlock_Lvals(6);
        %             if (lval05 == -1), par(i).refocpulse = 'lobw';
        %             else               par(i).refocpulse = 'exlobw'; end
        %             par(i).seqtype     = 'Voxel';
        %             par(i).TI          = h.WipMemBlock_Lvals(5);  % This is ALSO the Sat. Recovery time!!
        %             par(i).crusher     = h.WipMemBlock_Dvals(3);
        %             par(i).exciteppm   = h.WipMemBlock_Dvals(4);  % center of excite pulse (ppm)
        %             par(i).excitewidth = h.WipMemBlock_Dvals(5);  % width of excite pulse (ppm)
        %         case 3.0
        %             if (lval00 == -1), par(i).excitepulse = 'eburp';
        %             else               par(i).excitepulse = 'hsinc'; end
        %             lval05 = h.WipMemBlock_Lvals(6);
        %             if (lval05 == -1), par(i).refocpulse = 'lobw';
        %             else               par(i).refocpulse = 'exlobw'; end
        %             par(i).seqtype     = 'Voxel';
        %             par(i).TI          = h.WipMemBlock_Lvals(5);  % This is ALSO the Sat. Recovery time!!
        %             par(i).crusher     = h.WipMemBlock_Dvals(3);
        %             par(i).exciteppm   = h.WipMemBlock_Dvals(4);  % center of excite pulse (ppm)
        %             par(i).excitewidth = h.WipMemBlock_Dvals(5);  % width of excite pulse (ppm)
        %             par(i).presamp     = 2 * h.WipMemBlock_Lvals(8);   % samples before echo
        %         case 10.10
        %             rfpulses = {'eburp','hsinc'};
        %             index    = checkstruct(h,'WipMemBlock_Lval_02',0) + 1;
        %             par(i).seqtype       = 'Voxel';
        %             par(i).slicepolarity = 1 - 2 * h.WipMemBlock_Lvals(10);
        %             par(i).crushpolarity = 1 - 2 * h.WipMemBlock_Lvals(11);
        %             par(i).crushcycle    = double(checkstruct(h,'WipMemBlock_Lval_11',0));
        %             par(i).crusher       = h.WipMemBlock_Dvals(3);
        %             par(i).exciteppm     = h.WipMemBlock_Dvals(4);  % center of excite pulse (ppm)
        %             par(i).excitewidth   = h.WipMemBlock_Dvals(5);  % width of excite pulse (ppm)
        %             par(i).TI            = h.WipMemBlock_Lvals(7);     % This is ALSO the Sat. Recovery time!!
        %             par(i).presamp       = 2 * h.WipMemBlock_Lvals(8);   % samples before echo
        %             % handle TI/TS if no IR or SR done
        %             exptype = h.WipMemBlock_Lvals(6);     % experiment type, e.g. IR, SR, ...
        %             if (exptype == -1)
        %                 par(i).TI = 5000;  % essentially an M0 scan
        %             end
        %         case {10.11, 10.12, 10.13}
        %             lobes1   = h.WipMemBlock_Lvals(10);
        %             lobes2   = h.WipMemBlock_Lvals(11);
        %             rfpulses = {'eburp','hsinc',sprintf('UserSinc%1dx%1d',lobes1,lobes2)};
        %             index    = checkstruct(h,'WipMemBlock_Lval_02',0) + 1;
        %             par(i).excitepulse   = rfpulses{index};
        %             seqtypes = {'Voxel','Slice'};
        %             index    = checkstruct(h,'WipMemBlock_Lval_00',0) + 1;
        %             par(i).seqtype       = seqtypes{index};
        %             par(i).slicepolarity = 1 - 2 * checkstruct(h,'WipMemBlock_Lval_13',0);
        %             par(i).crushpolarity = 1 - 2 * checkstruct(h,'WipMemBlock_Lval_14',0);
        %             par(i).crushcycle    = double(checkstruct(h,'WipMemBlock_Lval_15',0));
        %             par(i).refocdur      = h.WipMemBlock_Lvals(5);
        %             par(i).crusher       = h.WipMemBlock_Dvals(3);
        %             par(i).exciteppm     = h.WipMemBlock_Dvals(4);  % center of excite pulse (ppm)
        %             par(i).excitewidth   = h.WipMemBlock_Dvals(5);  % width of excite pulse (ppm)
        %             par(i).TI            = h.WipMemBlock_Lvals(7);     % This is ALSO the Sat. Recovery time!!
        %             par(i).presamp       = 2 * h.WipMemBlock_Lvals(8);   % samples before echo
        %
        %             % Find experiment type
        %             explist  = {'svs', 'inv_selsinc', 'inv_selHS', 'inv_nonsel', 'sat_sel', 'sat_nonselAHP', 'sat_nonselBIR4', 'sat_nonselRECT', 'sat_nonselSINC'};
        %             expindex = double(checkstruct(h,'WipMemBlock_Lval_05',0)) + 1;
        %             par(i).exptype = explist{expindex};
        %             switch(par(i).exptype)
        %                 case 'svs'
        %                     par(i).TI = opt.m0_duration;  % essentially an M0 scan
        %             end
        %         otherwise
        %             fprintf(2,'ERROR: Unrecognized sequence version = %1.3f\n',par(i).seqver);
        %             return
        %     end
        % end
        %
        % % --- Phase cylce steps ---
        % switch par(i).seqtype
        %     case 'Voxel'
        %         par(i).phcyclesteps = 16;
        %     case 'Slice'
        %         par(i).phcyclesteps = 4;
        % end
    end

    par(i).status = 1;  % success
end
end

% ------------------------------------------------------------------------------------------------
% Get params from raw TWIX file
% ------------------------------------------------------------------------------------------------
function [par,npars] = raw_acqparams(datfiles, opt)

% -- make structure template ---
ntwix = numel(datfiles);
par   = make_par_struct(ntwix);

% --- now fill struct array ---
for i=1:ntwix
    % open file if given a filename
    if (ischar(datfiles{i}))
        if (opt.verbose), fprintf(1,'Reading raw TWIX file: %s\n',datfiles{i}); end
        [twix,rawfp] = me_mapVBVD(datfiles{i},'verbose',false); fclose(rawfp);
        if (iscell(twix)), twix = twix{2}; end
    else
        twix = datfiles{i};
    end

    % Handle multiple FIDs per rawdata file
    if (i == 1), nmeas = twix.image.NEco; end   % EJA returns multiple FIDs with this counter
    if (twix.image.NEco ~= nmeas), fprintf(2,'ERROR: Found multiple FIDs per TWIX file, but not the same number in each rawdata file (%1d and %1d)\n',nmeas,twix.image.NEco); return; end

    par(i).is_DFS       = contains(lower(twix.hdr.Phoenix.ProtocolName),'svssel') || contains(lower(twix.hdr.Phoenix.ProtocolName),'s3sel') || contains(lower(twix.hdr.Phoenix.ProtocolName),'selspec');
    par(i).is_EJA       = checkstruct(twix.hdr.Config,'Reps_eja',0);  % if this field exists, scan is from CMRR EJA
    par(i).is_CU        = contains(checkstruct(twix.hdr.Config,'SequenceString','xxx'),'_cu');

    par(i).nt           = checkstruct(twix.hdr.Meas,'RawCol',twix.hdr.Meas.VectorSize) * 2; % !!! NOTE OverSamping ON assumed !!!

    if (par(i).is_EJA)
        par(i).presamp  = twix.noise.iceParam(5); 
        %%par(i).presamp  = par(i).presamp * (~isempty(twix.hdr.Spice.RemoveOversampling) + 1);   % if the Dicoms will have OS removed, then EJA cuts the raw presamp val in half??!!
        par(i).postsamp = twix.image.NCol - par(i).presamp - par(i).nt;     % extra points at end?
    elseif (par(i).is_CU)
        par(i).presamp = twix.image.NCol - par(i).nt;
    end
    par(i).TR           = twix.hdr.Meas.TR/1E3;
    par(i).TE           = checkstruct2(twix.hdr.Meas,'TE','alTE',-1) / 1E3;  if (numel(par(i).TE) > 1), par(i).TE = par(i).TE(1); end
    par(i).f0           = twix.hdr.Meas.lFrequency / 1E6;
    par(i).dwelltime    = twix.hdr.Meas.alDwellTime(1) * 1E-9 * (opt.remove_oversamp + 1);  
    par(i).BW           = 1/par(i).dwelltime;
    parts               = split(twix.hdr.Meas.ExamMemoryUID,'_');
    par(i).scandate     = str2num(parts{4});
    par(i).scantime     = str2num(parts{5});  % Same for all series
    par(i).scanID       = twix.hdr.Phoenix.lSequenceID;
    par(i).nex          = checkstruct(twix.hdr.Config,'NumberOfAverages',twix.hdr.Config.Averages);
    par(i).actual_nex   = par(i).nex;
    par(i).coilname     = twix.hdr.Meas.TransmittingCoil;
    par(i).seriesname   = twix.hdr.Phoenix.ProtocolName;
    par(i).nucleus      = twix.hdr.Config.Nucleus;
    par(i).fieldstrength= [num2str(round(twix.hdr.Meas.flMagneticFieldStrength)) 'T'];
    par(i).txrefvolt    = twix.hdr.Meas.TransmitterReferenceAmplitude;
    par(i).voxdims      = [twix.hdr.Meas.VoiThickness, twix.hdr.Meas.VoiPhaseFOV, twix.hdr.Meas.VoiReadoutFOV];
    par(i).voxsize      = prod(par(i).voxdims)/1000;                     % voxel size (cc)
    par(i).seqbinary    = twix.hdr.Config.SequenceString;

    % --- Figure out DFS sequence version and wip params ---
    if (par(i).is_DFS)
        wipmem_struct = twix.hdr.Phoenix.sWipMemBlock.wipmem_struct;
        wipmem_long   = twix.hdr.Phoenix.sWipMemBlock.alFree;
        wipmem_dbl    = twix.hdr.Phoenix.sWipMemBlock.adFree;
        par(i)        = get_dfs_params(par(i),wipmem_struct,wipmem_long,wipmem_dbl,opt);
        par(i).postsamp = twix.image.NCol - par(i).presamp - par(i).nt;     % ONLY for raw data! - extra points at end?
    end
end
par(i).status = 1;  % success

% replicate par struct for each FID
if (nmeas > 1)
    par = repmat(par,nmeas,1); % par is now [nmeas,ntwix]
    for i=1:nmeas
        par(i,:).measnum = i;
    end
end
end


% ------------------------------------------------------------------------------------------------
% Figure out wipmem params from UPenn custom DFS sequence
% ------------------------------------------------------------------------------------------------
function par = get_dfs_params(par,wipmem_struct,wipmem_long,wipmem_dbl,opt)

osfactor = ~opt.remove_oversamp + 1;                % if remove OS, multiply or divide by 2

% --- Did we encode param names in the unused LONG values? ---
if (~isempty(wipmem_struct))
    seqtypes = {'Voxel','Slice'};
    rfpulses = {'eburp','hsinc',sprintf('UserSinc%1dx%1d',wipmem_struct.SncL,wipmem_struct.SncR)};
    explist  = {'svs', 'inv_selsinc', 'inv_selHS', 'inv_nonsel', 'sat_sel', 'sat_nonselAHP', 'sat_nonselBIR4', 'sat_nonselRECT', 'sat_nonselSINC'};
    par.seqver       = wipmem_struct.Vers;
    par.presamp      = wipmem_struct.DumP;
    par.excitepulse  = rfpulses{wipmem_struct.ExPu+1};
    par.seqtype      = seqtypes{wipmem_struct.LocM+1};
    par.exptype      = explist{wipmem_struct.InvM+1};
    par.presamp      = wipmem_struct.DumP * osfactor;
    par.adcdelay     = wipmem_struct.ADCD;
    par.exciteppm    = wipmem_struct.ExCn;
    par.excitewidth  = wipmem_struct.ExWd;
    par.TI           = wipmem_struct.TSTI;  if (wipmem_struct.InvM ==  0), par.TI = opt.m0_duration; end
    par.crushpolarity= wipmem_struct.RvCr * 2 - 1;
    par.slicepolarity= wipmem_struct.RvSS * 2 - 1;
    par.crushcycle   = wipmem_struct.CyCr;

% --- Otherwise, Complicated guesswork!! --
else
    wipmem_long = cellfun(@double, wipmem_long, 'UniformOutput', false); % convert longs to double
    wipmem_long = cell2mat(wipmem_long);
    wipmem_dbl  = cell2mat(wipmem_dbl);
    wipmem_long(isnan(wipmem_long)) = 0; %replacing NaNs w/ 0's
    wipmem_dbl(isnan(wipmem_dbl))   = 0; 

    if (wipmem_dbl(1) > 0 && wipmem_dbl(1) < 1)                % this means there is an explicit version variable
        par.seqver = 10 + wipmem_dbl(1);
    else
        fprintf(2,'WARNING: Cant figure DFS sequence version from wipmem_dbl(1) = %g and wipmem_long(1) = %1d\n',wipmem_dbl(1),wipmem_long(1));
        return
    end

    % --- Assign vars from WipMem ---
    switch(par.seqver)
        case 10.10
            par.seqtype     = 'Voxel';
            par.presamp = wipmem_long(8) * osfactor;
        case {10.11, 10.12, 10.13}
            par.presamp = wipmem_long(8) * osfactor;
            par.adcdelay = wipmem_long(9);
            par.exciteppm     = wipmem_dbl(4);
            par.excitewidth   = wipmem_dbl(5);
            par.TI            = wipmem_long(7);
            par.crushcycle    = wipmem_long(16);
            explist  = {'svs', 'inv_selsinc', 'inv_selHS', 'inv_nonsel', 'sat_sel', 'sat_nonselAHP', 'sat_nonselBIR4', 'sat_nonselRECT', 'sat_nonselSINC'};
            expindex = wipmem_long(6)  + 1;
            par.exptype = explist{expindex};
            switch(par.exptype)
                case 'svs'
                    par.TI = opt.m0_duration;  % essentially an M0 scan
            end
            lobes1   = wipmem_long(10);
            lobes2   = wipmem_long(11);
            rfpulses = {'eburp','hsinc',sprintf('UserSinc%1dx%1d',lobes1,lobes2)};
            index    = wipmem_long(3) + 1;
            par.excitepulse   = rfpulses{index};

            seqtypes = {'Voxel','Slice'};
            index    = wipmem_long(1) + 1;
            par.seqtype       = seqtypes{index};
        otherwise
            fprintf(2,'ERROR: Unrecognized sequence version = %1.3f\n',par.seqver);
            return
    end
end

% --- Phase cylce steps ---
switch par.seqtype
    case 'Voxel'
        par.phcyclesteps = 16;
    case 'Slice'
        par.phcyclesteps = 4;
end
end



% ------------------------------------------------------------------------------------------------
% Make default param structure
% ------------------------------------------------------------------------------------------------
function par = make_par_struct(nstruct)

if (nargin < 1), nstruct = 1; end
par.status          = -1;
par.TR              = -1;
par.TE              = -1;
par.nt              = -1;
par.nex             = -1;
par.actual_nex      = -1;
par.TI              = -1;
par.crusher         = -1;
par.exciteppm       = -99;
par.excitewidth     = -1;
par.dwelltime       = -1;
par.BW              = -1;
par.f0              = -1;
par.txrefvolt       = -1;
par.seriesnum       = -1;
par.imagenum        = -1;
par.measnum         = -1;  % for rawdata files w/ multiple FIDs, mainly EJA
par.seriesname      = 'unknown';
par.voxdims         = [-1, -1, -1];
par.voxsize         = -1;
par.scandate        = -1;
par.scantime        = -1;
par.patname         = 'anon';
par.patage          = 0;
par.coilname        = '';
par.excitepulse     = 'unknown';
par.refocdur        = -1;
par.refocpulse      = '';
par.seqtype         = 'svs';
par.crushpolarity   = -99;
par.slicepolarity   = -99;
par.crushcycle      = 0;
par.presamp         = 0;
par.postsamp        = 0;
par.adcdelay        = 0;
par.phcyclesteps    = 0;
par.seqver          = -1;
par.scanID          = -1;
par.is_DFS          = 0;    % CAMIPM svssel
par.is_EJA          = 0;    % CMRR EJA
par.is_CU           = 0;    % Columbia
par.exptype         = 'unknown';
par.seqbinary       = 'unknown';
if (nstruct > 1), par = repmat(par,1,nstruct); end
end




