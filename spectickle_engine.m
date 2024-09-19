% ----------------------------------------------------------------------------------
% --- SPECTICKLE_ENGINE() ---
% ----------------------------------------------------------------------------------
function [status,res,opt,par,data] = spectickle_engine(dfsfiles,opt_dfs,waterfiles,opt_water,peaktable,exampdicom,comps)

status = 0; % means error 
res  = [];
opt  = [];
par  = [];
data = [];

% --- Check calling syntax ---
if (nargin < 2) || (nargin == 3)
    fprintf(2,'Calling syntax error: %s(svsfiles,opt_svs,[reffiles,opt_ref],[example_dicomfile],[comps])\n',mfilename());
    return
end
if (nargin < 3), waterfiles = {}; end
if (nargin < 5), peaktable  = []; end
if (nargin < 6), exampdicom = []; end
if (nargin < 7), comps      = []; end

st = dbstack; subroutine_name = st.name;
fprintf(1,'\n-------------------------------------------------------------\n');
fprintf(1,'-------------------------------------------------------------\n');
fprintf(1,'Beginning %s()\n',subroutine_name);
fprintf(1,'-------------------------------------------------------------\n');
fprintf(1,'-------------------------------------------------------------\n');

% --- Check options ---
opt_dfs   = check_options(opt_dfs);
opt_water = check_options(opt_water);
if (~iscell(dfsfiles)),   dfsfiles   = {dfsfiles}; end
if (isempty(waterfiles)), waterfiles = {}; end
if (~iscell(waterfiles)), waterfiles = {waterfiles}; end

% --- Get data folder ---
tmp    = fileparts(dfsfiles{1});
topdir = fileparts(tmp); 
[~,opt_dfs.data_name] = fileparts(topdir); % get folder name to use as dataset identifier

% --- Load peaktable (if not provided as a structure, then it can be the filename, or empty string) ---
if (~isstruct(peaktable))
    [peaktable,peakfile] = spectickle_read_peaktable(peaktable); 
    if (isempty(peaktable)), fprintf(2,'ERROR: Could not find/load peak table file "%s"\n',peakfile); return; end
end

% --- Build peaks structures for peaks to fit ---
peaks_dfs = spectickle_get_peaks(peaktable,opt_dfs); if (isempty(peaks_dfs)), return; end
if (~isempty(waterfiles))
    peaks_water = spectickle_get_peaks(peaktable,opt_water); if (isempty(peaks_water)), return; end
end

% ----------------------------------
% --- Read Dicom files ---
% ----------------------------------
ndfs = numel(dfsfiles);
nwat = numel(waterfiles);
if (~opt_dfs.use_rawdata) 
    % --- Read Dicom water reference files ---
    if (~isempty(waterfiles))
        hdr_water = dicom_header(waterfiles);
        par_water = svs_acqparams(hdr_water,opt_water);
        for i=1:nwat
            fid = siemens_dicom_mrs(hdr_water{i},opt_water);
            if (i==1), fids_water = repmat(fid,1,nwat); end
            fids_water(:,i) = fid;
        end
    end
    % --- Read Dicom "DFS" files ---
    hdr_dfs = dicom_header(dfsfiles);
    par_dfs = svs_acqparams(hdr_dfs,opt_dfs);
    for i=1:ndfs
        if (i==1)
            [fid,h,t,hz,ppm] = siemens_dicom_mrs(hdr_dfs{i},opt_dfs);
            fids_dfs = repmat(fid,1,ndfs);
        else
            fid = siemens_dicom_mrs(hdr_dfs{i},opt_dfs);
        end
        fids_dfs(:,i) = fid;
    end
    combine_source = 'dicom';    % remember that coil combine was done by ICE via Dicoms

% -----------------------------------
% --- Read TWIX rawdata files ---
% -----------------------------------
else
    opt_dfs.raw_scale = 1E6;
    opt_water.raw_scale = opt_dfs.raw_scale;
    % --- Read TWIX water reference files ---
    if (~isempty(waterfiles))
        [fids_water,~,coeffs,par_water, ~,~,~,fshift,whitemat,chanphase] = siemens_raw_mrs(waterfiles,[],[],[],[],opt_water);
        if (isempty(fids_water)); return; end
        combine_source = 'water';    % remember that coil combine was from water reference
    else
        coeffs   = [];
        fshift   = [];
        whitemat = [];
        chanphase = [];
        combine_source = 'self';    % remember that coil combine was from DFS files        
    end
    % --- Read TWIX "DFS" files ---
    [fids_dfs,~,~,par_dfs,t,hz,ppm] = siemens_raw_mrs(dfsfiles,coeffs,fshift,whitemat,chanphase,opt_dfs);
    if (isempty(fids_dfs)); return; end
    [np,nmeas,nfiles] = size(fids_dfs);
    fids_dfs = reshape(fids_dfs,np,nmeas*nfiles); 
    par_dfs  = reshape(par_dfs,nmeas*nfiles,1);
    
    % --- Fill in some params which might not exist in raw TWIX files ---
    if (~isempty(exampdicom))
        par_dcm = svs_acqparams(exampdicom,opt_dfs);
        if (par_dfs(1).scanID ~= par_dcm.scanID)
            fprintf(2,'ERROR: Raw data scanID (%1d) != Dicom scanID (%1d)\n',par_dfs(1).scanID,par_dcm.scanID); 
            % return;
        end
        if (par_dfs(1).scandate ~= par_dcm.scandate) || (abs(par_dfs(1).scantime-par_dcm.scantime) > (60*60*3))
            fprintf(2,'ERROR: Raw data scandate/scantime (%1d/%1d) != Dicom scandate/scantime (%1d/%1d)\n',par_dfs(1).scandate,par_dfs(1).scantime,par_dcm.scandate,par_dcm.scantime); 
            return;
        end
        for i=1:ndfs, par_dfs(i).patname = par_dcm.patname; end
        for i=1:ndfs, par_dfs(i).patage  = par_dcm.patage; end
    end
end

% Combine subtraction pairs for spectral editing
if (opt_dfs.subtract_pairs)    
    [np,ndfs_fids] = size(fids_dfs);
    if (ndfs_fids == 1),        fprintf(2,'ERROR: Option for spectral editing is set, but only 1 FID provided.\n'); return; end
    if (mod(ndfs_fids,1) ~= 0), fprintf(2,'ERROR: Option for spectral editing is set, but an odd number (%1d) of FIDs provided.\n',ndfs_fids); return; end
    fprintf(1,'Spectral editing set. Subtracting %1d consecutive pairs of FIDs.\n',ndfs_fids/2);
    fids_dfs = reshape(fids_dfs,np,ndfs_fids/2,2);
    fids_dfs = fids_dfs(:,:,2) - fids_dfs(:,:,1); % subtract pairs of FIDs
    if (ndfs ~= (ndfs_fids/2))  % e.g. 2 Dicom files per editing pair, decimate info structs by 2x
        fprintf(1,'\tDecimating info structs by 2x\n');
        index    = 1:2:(ndfs);
        par_dfs  = par_dfs(index);
        dfsfiles = dfsfiles(index);
    end
end

% handle different number of FIDs than data files (e.g. multiple measurements/file)
ndfs_fids = size(fids_dfs,2);  % get this again in case spectral editing was done
if (ndfs_fids ~= ndfs)
    nmeas = ndfs_fids/ndfs;  % this is how many FIDs were in each data file
    for i=1:ndfs
        for j=1:nmeas
            new_dfsfiles{j,i} = dfsfiles{i};  % replicate the filenames so that each FID is associated with a filename
        end
    end
    dfsfiles = reshape(new_dfsfiles,nmeas*ndfs,1);
end
ndfs = ndfs_fids; % this is the correct number of FIDs now

% ----------------------------------------------------------------------------------
% --- PROCESS WATER FILES ---
% ----------------------------------------------------------------------------------
if (~isempty(waterfiles))
    if (opt_water.verbose), fprintf(1,'Processing %1d Water spectra...\n',nwat); end
    opt_water.data_name      = opt_dfs.data_name;
    opt_water.baseline_lwmin = 100;
    
    % --- Sort ---
    if (nwat > 1)
        opt_water.file_labels = cell(1,nwat);  % labels for each data set (for plot titles)
        switch lower(opt_water.sort_by)
            case 'series'
                fprintf(1,'Sorting Water files by series number\n');
                [~,order] = sort([par_water.seriesnum]);
                for i=1:nwat, opt_water.file_labels{i} = sprintf('Series %1d',par_water(order(i)).seriesnum); end
            case 'te'
                fprintf(1,'Sorting Water files by TE\n');
                [~,order] = sort([par_water.TE]);
                for i=1:nwat, opt_water.file_labels{i} = sprintf('TE = %1.1f',par_water(order(i)).TE); end
            case 'tr'
                fprintf(1,'Sorting Water files by TR\n');
                [~,order] = sort([par_water.TR],'descend');
                for i=1:nwat, opt_water.file_labels{i} = sprintf('TR = %1.0f',par_water(order(i)).TR); end
            case 'ti'
                fprintf(1,'Sorting Water files by TI\n');
                [~,order] = sort([par_water.TI],'descend');
                for i=1:nwat, opt_water.file_labels{i} = sprintf('TI = %1.0f',par_water(order(i)).TI); end
            case {'sr','ts'}
                fprintf(1,'Sorting Water files by TS\n');
                [~,order] = sort([par_water.TI],'descend');
                for i=1:nwat, opt_water.file_labels{i} = sprintf('TS = %1.0f',par_water(order(i)).TI); end
            case 'nex'
                fprintf(1,'Sorting Water files by number of averages\n');
                [~,order] = sort([par_water.nex],'descend');
                for i=1:nwat, opt_water.file_labels{i} = sprintf('Navg = %1d',par_water(order(i)).nex); end
            otherwise
                fprintf(1,'No sorting of Water files done\n');
                order = 1:nwat;
                for i=1:nwat, opt_water.file_labels{i} = sprintf('Scan %1d',i); end
        end
        par_water  = par_water(order);
        waterfiles = waterfiles(order);
        fids_water = fids_water(:,order); % !!!!
    else
        opt_water.file_labels = {sprintf('Series %1d',par_water.seriesnum)};
    end
    
    % --- Do HSVD fit ---
    if (opt_water.separate_model)
        opt_water.filenames = waterfiles;  % this is just so the hsvd_engine() can print what data file it is processing
        for i=1:nwat
            [res_water(i),opt_water,data_water] = hsvd_engine(fids_water(:,i),opt_water,par_water(i),peaks_water,t,hz,ppm);
            res_water(i).combine_source = combine_source;
            if (res_water(i).status == 0), return; end
        end
    else
        [res_water,opt_water,data_water] = hsvd_engine(fids_water,opt_water,par_water,peaks_water,t,hz,ppm);
        res_water(:).combine_source = combine_source;
        if (~all([res_water.status])), return; end
    end
    opt_dfs.use_ppm = opt_water.ppm_used;   % if needed, use ppm set from water ref
else
    res_water.status = -1;              % for data save later, indicate no water ref peak found
    opt_water.status = -1;
    par_water.status = -1;
    data_water.status = -1;
end

% ----------------------------------------------------------------------------------
% --- PROCESS DFS SPECTRA ---
% ----------------------------------------------------------------------------------
if (opt_dfs.verbose), fprintf(1,'Processing %1d SVS spectra...\n',ndfs); end
% --- Sort ---
if (ndfs > 1)
    opt_dfs.file_labels = cell(1,ndfs);  % labels for each data set (for plot titles)
    switch lower(opt_dfs.sort_by)
        case 'series'
            fprintf(1,'Sorting SVS files by series number\n');
            [~,order] = sort([par_dfs.seriesnum]);
            for i=1:ndfs, opt_dfs.file_labels{i} = sprintf('Series %1d: %s',par_dfs(order(i)).seriesnum,par_dfs(order(i)).seriesname); end
        case 'te'
            fprintf(1,'Sorting SVS files by TE\n');
            [~,order] = sort([par_dfs.TE]);
            for i=1:ndfs, opt_dfs.file_labels{i} = sprintf('TE = %1.1f',par_dfs(order(i)).TE); end
        case 'tr'
            fprintf(1,'Sorting SVS files by TR\n');
            [~,order] = sort([par_dfs.TR],'descend');
            for i=1:ndfs, opt_dfs.file_labels{i} = sprintf('TR = %1.0f',par_dfs(order(i)).TR); end
        case 'ti'
            fprintf(1,'Sorting SVS files by TI\n');
            [~,order] = sort([par_dfs.TI],'descend');
            for i=1:ndfs, opt_dfs.file_labels{i} = sprintf('TI = %1.0f',par_dfs(order(i)).TI); end
        case {'sr','ts'}
            fprintf(1,'Sorting SVS files by TS\n');
            [~,order] = sort([par_dfs.TI],'descend');
            for i=1:ndfs, opt_dfs.file_labels{i} = sprintf('TS = %1.0f',par_dfs(order(i)).TI); end
        case 'nex'
            fprintf(1,'Sorting SVS files by number of averages\n');
            [~,order] = sort([par_dfs.nex],'descend');
            for i=1:ndfs, opt_dfs.file_labels{i} = sprintf('Navg = %1d',par_dfs(order(i)).nex); end
        case 'measurement'
            fprintf(1,'Sorting SVS files by measurement number\n');
            [~,order] = sort([par_dfs.imagenum]);
            for i=1:ndfs, opt_dfs.file_labels{i} = sprintf('Measurement = %1d',par_dfs(order(i)).imagenum); end
        otherwise
            fprintf(1,'No sorting of SVS files done\n');
            order = 1:ndfs;
            for i=1:ndfs, opt_dfs.file_labels{i} = sprintf('Scan %1d',i); end
    end
    par_dfs  = par_dfs(order);
    dfsfiles = dfsfiles(order);
    fids_dfs = fids_dfs(:,order); % !!!!
else
%    opt_dfs.file_labels = {sprintf('%s (series #%1d)',par_dfs.seriesname,par_dfs.seriesnum)};
    opt_dfs.file_labels = {par_dfs.seriesname};
end

% --- Freq drift correct the FIDs ---
if (ndfs > 1) && (opt_dfs.fid_freqalign)
    if (opt_dfs.verbose), fprintf(1,'Performing frequency drift correction between scans (reference scan is "%s")\n',opt_dfs.file_labels{1}); end
    fidref = fids_dfs(:,1); 
    np     = size(fidref,1);
    for i=2:ndfs
        [fshift,ffit,fpars] = fidfreqshift(t,fids_dfs(:,i),fidref,[0,1,0],np/2);
        if (opt_dfs.verbose), fprintf(1,'Frequency drift correction of %4.1f Hz applied to scan "%s"\n',fshift,opt_dfs.file_labels{i}); end
        fids_dfs(:,i) = ffit / fpars(2);  % undo amplitude fit
    end
end

% --- Separately Model and Fit blocks of FIDS ---
if (opt_dfs.separate_model > 0)
    nblocks = ceil(ndfs/opt_dfs.separate_model);
    for i=1:nblocks
        index = (1:opt_dfs.separate_model) + (i-1)*opt_dfs.separate_model;
        index = index(index <= ndfs);
        opt_dfs.model_avgfids = [1 1000]; % this mean sums whole block for Modeling step
        opt_dfs.filenames = dfsfiles(index);
       
        [res_dfs(i),opt_dfsX(i),data_dfs(i)] = hsvd_engine(fids_dfs(:,index),opt_dfs,par_dfs(index),peaks_dfs,t,hz,ppm);
        if (res_dfs(i).status == 0), return; end
        res_dfs(i).combine_source = combine_source;
    end
    opt_dfs = opt_dfsX; % problem w/ changing input opt_dfs during loop above

% --- Model and Fit all FIDS together ---
else
    opt_dfs.filenames = dfsfiles;
    [res_dfs,opt_dfs,data_dfs] = hsvd_engine(fids_dfs,opt_dfs,par_dfs,peaks_dfs,t,hz,ppm,comps);
    for j=1:numel(res_dfs), res_dfs(j).combine_source = combine_source; end
    if (~all([res_dfs.status])), return; end
end

% --- Bundle up results ---
res.res_dfs1   = res_dfs;  res.res_dfs2   = []; res.res_water   = res_water;
opt.opt_dfs1   = opt_dfs;  opt.opt_dfs2   = []; opt.opt_water   = opt_water;
par.par_dfs1   = par_dfs;  par.par_dfs2   = []; par.par_water   = par_water;
data.data_dfs1 = data_dfs; data.data_dfs2 = []; data.data_water = data_water;
data.peaktable = peaktable;

status = 1; % success

fprintf(1,'-------------------------------------------------------------\n');
fprintf(1,'Done %s()\n',subroutine_name);
fprintf(1,'-------------------------------------------------------------\n\n');
end

% -------------------------------------------------------------------------------------
% --- HSVD_ENGINE() --
% -------------------------------------------------------------------------------------
function [res,opt,data] = hsvd_engine(fids,opt,par,peaks,t,hz,ppm,comps)
res = [];
% --- handle inputs/outputs ---
savedata = 0;
if (nargin  < 8), comps = []; end
if (nargout > 2)
    savedata    = 1; 
    data.status = -1;
end
bw = 1/(t(2)-t(1));

st = dbstack(); subroutine_name = st.name;
fprintf(1,'\n-------------------------------------------------------------\n');
fprintf(1,'Beginning %s()\n',subroutine_name);
fprintf(1,'-------------------------------------------------------------\n');
if (opt.verbose), fprintf(1,'Processing study "%s" for experiment type "%s"\n',opt.data_name,opt.exptype); end

% --------------------------------------------------------------------------------------
% --- DEBUGGING CODE - Truncate FID readout to equal HSVD model points
% --------------------------------------------------------------------------------------
if (opt.debug && opt.hsvd_np > 0)
    fprintf(2,'Truncating data to have only %1d FID points\n',opt.hsvd_np)
    fids       = fids(1:opt.hsvd_np,:);
    [t,hz,ppm] = spec_xaxis(opt.hsvd_np,bw,par.f0);
end

% --- initialize returned items ---
[np,nfids] = size(fids);
switch lower(opt.avg_group)  % handle averaging by group
    case 'none'
        nout = nfids;
    case 'all'
        nout = 1;
    case 'tr'
        avgvals = unique([par.TR]);
        nout    = numel(avgvals);
    case 'te'
        avgvals = unique([par.TE]);
        nout    = numel(avgvals);
    case {'ti','ts'}
        avgvals = unique([par.TI]);
        nout    = numel(avgvals);
    case 'series'
        avgvals = unique([par.seriesnum]);
        nout    = numel(avgvals);
    case {'every','every...'}
        nout = ceil(nfids/opt.fit_avgfids);
    otherwise
        fprintf(2,'ERROR: Dont know how to average FIDs by %s\n',opt.avg_group); return; 
end
res = create_res_struct(nout);

% --- deal with peaks to find ---
peaknames = [peaks.labels;'baseline';'unassigned'];    % append 'baseline' or 'unassigned' for all components not assigned to a peak
ppm0 = min(peaks.ppm) - 1.;  % ppm of interest range
ppm1 = max(peaks.ppm) + 1.;

% --- If labels for dicom files not specified... ---
if (isempty(opt.file_labels))
    opt.file_labels = cell(1,nfids);
    for i=1:nfids
        opt.file_labels{i} = '';
    end
end

% --- Get FID(s) to be MODELED by HSVD ---
if (opt.model_avgfids(2) > nfids), opt.model_avgfids(2) = min(opt.model_avgfids(2),nfids); end  % make 1:nfids at most
if (opt.model_avgfids(1) > opt.model_avgfids(2)), fprintf(2,'ERROR: model_avgfids has element2 > element1!\n'); return; end
navg = opt.model_avgfids(2)-opt.model_avgfids(1)+1;
if (navg > 1 && opt.verbose), fprintf(1,'Averaging %1d FIDs (from %1d to %1d)...\n',navg,opt.model_avgfids(1),opt.model_avgfids(2)); end
if (opt.verbose)
    fprintf(1,'Modeling FIDs:\n')
    for i=opt.model_avgfids(1):opt.model_avgfids(2)
        [path,fname,ext] = fileparts(opt.filenames{i});
        [~,tpath] = fileparts(path);
        fprintf(1,'    %s\n',['..' filesep() tpath filesep() fname ext]);
    end
end
fid = mean(fids(:,opt.model_avgfids(1):opt.model_avgfids(2)),2);

% --- Decide how to center the ppm scale ---
switch lower(opt.center_ppm)
    case 'center'
        if (opt.verbose), fprintf(1,'Setting ppm axis to %1.1f ppm at spectrum center\n',opt.ppm_center); end
        ppm      = ppm + opt.ppm_center;     
    case 'max'
        if (opt.verbose), fprintf(1,'Setting ppm axis to %1.1f ppm at spectrum max\n',opt.ppm_center); end
        spec     = m1dfft(fid);
        [~,imax] = max(abs(spec(:)));
        ppmshift = ppm(imax);
        ppm      = ppm - ppmshift + opt.ppm_center;
    case {'external','water'}
        if (opt.verbose), fprintf(1,'Using provided ppm scale\n'); end
        ppm = opt.use_ppm;
    otherwise
        fprintf('ERROR: Unrecognized center_ppm option: %s\n',opt.center_ppm);
        return
end
ppmshift = ppm(np/2+1); % final center PPM we used
if (opt.verbose), fprintf(1,'Center of ppm axis is %1.1f\n',ppmshift); end
opt.ppm_used = ppm;     % save what ppm scale was

% --- Line broaden  ---
if (opt.lb > 0.)
    if (opt.verbose), fprintf(1,'Applying %1.1f Hz line broadening\n',opt.lb); end
    apod = exp(-t*opt.lb*pi);
    fid  = fid.*apod;
end

% --------------
% --- HSVD!! ---
% --------------
if (isempty(comps))
    if (opt.hsvd_np == 0)   % default choice for HSVD number of model points
        if (np >= 1024), opt.hsvd_np = 1024;
        else,            opt.hsvd_np = np;   end  
    end                    % user chose the HSVD model points size
    if (opt.hsvd_np < (2*opt.hsvd_order)), fprintf(2,'ERROR: Number of FID points to model (%1d) must be > 2 * HSVD order (%1d)\n',opt.hsvd_np,opt.hsvd_order); return; end
    if ~mod(opt.hsvd_np,2), opt.hsvd_np = opt.hsvd_np - 1; end    % odd number makes the Hankel matrix square (doesn't seem to matter)

    if (opt.verbose), fprintf(1,'Performing %s on %1d time points for %1d components\n',opt.model_algo,opt.hsvd_np,opt.hsvd_order); end
    switch opt.model_algo
        case 'HSVD'
            [fitfid,basis,pars,~,amp] = me_hsvd(fid,opt.hsvd_np,opt.hsvd_order,bw,opt.fit_skippts,opt.verbose,0);
        case 'HLSVD'
            [fitfid,basis,pars,~,amp] = me_hsvd(fid,opt.hsvd_np,opt.hsvd_order,bw,opt.fit_skippts,opt.verbose,1);
        case 'HTLS'
        [fitfid,basis,pars,~,amp] = me_htls(fid,opt.hsvd_np,opt.hsvd_order,bw,opt.fit_skippts,opt.verbose);
    end
    nbasis = size(basis,2);
    if (opt.verbose), fprintf(1,'\tFound %1d valid components\n',nbasis); end
    
    % --- calc LWs in Hz and freqs in ppm ---
    lw       = 1000./(pi*pars(:,3));   % Hz
    freq_ppm = -pars(:,1)/par(1).f0 + ppmshift;
    
    % --- display initial params ---
    if (opt.verbose)
        dump_comps(1:nbasis,freq_ppm,amp,lw,pars(:,4));
    end
    
    % ------------------------------------
    % --- weed out unwanted components ---
    % ------------------------------------
    amp_rel = abs(amp)/max(abs(amp(:)));
    pamp  = (amp_rel  >= opt.amp_min);                 nkeep_amp  = sum(pamp);   % this is a logical array, not an index
    plw1  = (lw       <= opt.lw_max);                  nkeep_lw1  = sum(plw1);
    plw2  = (lw       >= opt.lw_min);                  nkeep_lw2  = sum(plw2);
    pfreq = (freq_ppm >= opt.ppm_min);                 nkeep_freq = sum(pfreq);
    
    if (nkeep_amp  < nbasis) && (opt.verbose), fprintf(1,'Throwing out %1d components with relative amplitude < %1.4f\n',nbasis-nkeep_amp,opt.amp_min); end
    if (nkeep_lw1  < nbasis) && (opt.verbose), fprintf(1,'Throwing out %1d components with LW > %1.2f\n',nbasis-nkeep_lw1,opt.lw_max); end
    if (nkeep_lw2  < nbasis) && (opt.verbose), fprintf(1,'Throwing out %1d components with LW < %1.2f\n',nbasis-nkeep_lw2,opt.lw_min); end
    if (nkeep_freq < nbasis) && (opt.verbose), fprintf(1,'Throwing out %1d components with ppm freq < %1.2f\n',nbasis-nkeep_freq,opt.ppm_min); end
    
    p = find(pamp & plw1 & plw2 & pfreq); nkeep = numel(p);                    % this is an index
    basis    = basis(:,p);
    pars     = pars(p,:);
    amp      = amp(p);
    lw       = lw(p);
    freq_ppm = freq_ppm(p);
    nbasis   = nkeep;
    
    % --- prepare structure with final HSVD results and fits to the data ---
    comps.ncomps = nbasis;
    comps.basis  = complex(zeros(np,comps.ncomps),zeros(np,comps.ncomps));
    for i = 1:comps.ncomps
        comps.basis(:,i) = basis(:,i)/amp(i); % normalize basis amps AND zero their PHASE!!
    end
    comps.basis_inv = pinv(comps.basis((opt.fit_skippts+1):end,:));      % gonna need this a lot (NOTE skipping bad initial FID points)
    comps.freq_hz  = pars(:,1);               % Hz
    comps.freq_ppm = freq_ppm;
    comps.t2star   = pars(:,3);               % ms
    comps.phase    = pars(:,4);               % deg
    comps.lw       = lw;   % Hz
    
    % --- assign components to peaks ---
    comps.peak_assignment = zeros(comps.ncomps,1) + peaks.count + 1;    % default assignment is 'baseline'
    comps.missing_peaks   = zeros(peaks.count,1,'logical'); % will hold boolean list of missing peaks
    for i=1:peaks.count
        p = find((comps.lw  >= peaks.lw_min(i))                    & ...
            (comps.lw       <= peaks.lw_max(i))                    & ...
            (comps.phase    <= peaks.phase_max(i))                 & ...
            (comps.phase    >= peaks.phase_min(i))                 & ...
            (comps.freq_ppm >= (peaks.ppm(i) - peaks.ppm_tol(i)))  & ...
            (comps.freq_ppm <= (peaks.ppm(i) + peaks.ppm_tol(i))));
        if (~isempty(p))
            comps.peak_assignment(p) = i;
        else
            fprintf(2,'WARNING: Could not find any components for %s peak\n',peaks.labels{i});
            comps.missing_peaks(i) = true; % remember
        end
    end
    
    % --- Weed out narrow peaks from baseline and call them 'unassigned' ---
    if (opt.baseline_lwmin > 0)
        p = find((comps.lw < opt.baseline_lwmin) & (comps.peak_assignment == (peaks.count + 1)));
        if (numel(p) > 0), comps.peak_assignment(p) = peaks.count + 2; end
    end
    
    % --- For print-outs, just show components in ppm range of interest
    pdisp = find((comps.freq_ppm >= ppm0) & (comps.freq_ppm <= ppm1)); % only show components within spec range of interest
    ndisp = numel(pdisp);
    
    % --- print component list w/ peak assignments ---
    if (opt.verbose)
        fprintf(1,'%1d components are within ppm range of interest [%1.2f,%1.2f]\n',ndisp,ppm0,ppm1);
        dump_comps(pdisp,comps.freq_ppm(pdisp),amp(pdisp),comps.lw(pdisp),pars(pdisp,4),{peaknames{comps.peak_assignment(pdisp)}});
    end
    
    % --- Handle if any peaks were not assigned at least one component ---
    if (all(comps.missing_peaks)),                            fprintf(2,'ERROR: Could not find ANY peaks!\n'); return; end
    if (any(comps.missing_peaks)) && (~opt.missing_peaks_ok), fprintf(2,'ERROR: Could not find all peaks.\n'); return; end
    
    % --- Keep only one component per peak? ---
    if (~isempty(opt.single_component))
        if (ischar(opt.single_component)) % if just one entry, e.g. 'ppm' repeat it for all peaks
            tmp = opt.single_component;
            opt.single_component = cell(1,peaks.count);
            opt.single_component(:) = {tmp};
        end
        limit_implemented = 0;
        for i=1:peaks.count
            p = find(comps.peak_assignment == i);
            if (numel(p) > 1)
                switch lower(opt.single_component{i})
                    case {'','off','-- off --'}             % do nothing
                    case {'amp','amplitude'}
                        [~,order] = sort(amp(p));
                        comps.peak_assignment(p(order(1:end-1))) = peaks.count + 1;    % assign smaller components back to baseline
                        if (opt.verbose), fprintf(1,'Limiting %s peak to a single component w/ largest amplitude\n',peaks.labels{i}); end
                        limit_implemented = 1;
                    case {'ppm','frequency','freq'}
                        [~,order] = sort(abs(comps.freq_ppm(p) - peaks.ppm(i)));
                        comps.peak_assignment(p(order(2:end))) = peaks.count + 1;    % assign further away components back to baseline
                        if (opt.verbose), fprintf(1,'Limiting %s peak to a single component closest in ppm to peak template\n',peaks.labels{i}); end
                        limit_implemented = 1;
                    otherwise
                        fprintf(2,'ERROR: Unknown "single_component" option: "%s"\n',opt.single_component{i});
                        return
                end
            end
        end
        
        % --- Weed out narrow peaks (again) from baseline and call them 'unassigned' ---
        if (opt.baseline_lwmin > 0)
            p = find((comps.lw < opt.baseline_lwmin) & (comps.peak_assignment == (peaks.count + 1)));
            if (numel(p) > 0), comps.peak_assignment(p) = peaks.count + 2; end
        end
        
        if (opt.verbose && limit_implemented)
            fprintf(1,'Peak assigments after implementing "single_component" option\n');
            dump_comps(pdisp,comps.freq_ppm(pdisp),amp(pdisp),comps.lw(pdisp),pars(pdisp,4),{peaknames{comps.peak_assignment(pdisp)}});
        end
    end

    % if (opt.yoke_peakgroups)
    %     if (all(matches(peaks.yoke_group,'none')))
    %         fprintf(2,'WARNING: Yoking option selected but no yoked peak groups were selected\n');
    %     else
    %         p = ~matches(peaks.yoke_group,'none');  % ignore peaks in yoke group 'none'
    %         groups = unique(peaks.yoke_group(p));   % run through each yoke group
    %         for i=1:numel(groups)
    %             p = find(matches(peaks.groups,groups{i}));  % find peaks in this yoke group
    %             for j=1:numel(p)
    %                 index = find(matches(peaknames,peaks.labels(p(j))));  % this is the peak assignment value for this peak.label
    %                 pcomp = find(comps.peak_assignment == index);
    % 
    % 
    %                 xx=5;
    % 
    %             end
    % 
    %         end
    %     end
    % end

    % --- save data results ---
    if (savedata)
        data.modelfid    = fid;
        data.modelfidfit = fitfid;
    end

% --- Use HSVD basis set provided in call ---
else
    if (opt.verbose), fprintf(2,'NOTE: Re-using provided Model basis set with %1d components\n',comps.ncomps); end
    %pdisp = find((comps.freq_ppm >= ppm0) & (comps.freq_ppm <= ppm1)); % only show components within spec range of interest
end

% ------------------------------------------------------
% --- Regression to HSVD basis set  for amps/phases ---
% ------------------------------------------------------
comps.amp = complex(zeros(comps.ncomps,nout),zeros(comps.ncomps,nout));    
for i=1:nout
    % --- Get index of FIDs to average for this output block ---
    switch lower(opt.avg_group)
        case 'none'
            index = i;
        case 'all'
            if (opt.verbose), fprintf(1,'Averaging ALL %1d FIDs for regression to HSVD basis set\n',nfids); end
            index = 1:nfids;
            opt.file_labels{i} = sprintf('Average of Scans %1d:%1d',index(1),index(end)); % overwrite fit result label
        case 'tr'
            index = find([par.TR] == avgvals(i));
            opt.file_labels{i} = sprintf('TR = %1.1f',avgvals(i)); % overwrite fit result label
            if (opt.verbose), fprintf(1,'Averaging %1d FIDs with TR=%1.1f for regression to HSVD basis set\n',numel(index),avgvals(i)); end
        case 'te'
            index = find([par.TE] == avgvals(i));
            opt.file_labels{i} = sprintf('TE = %1.1f',avgvals(i)); % overwrite fit result label
            if (opt.verbose), fprintf(1,'Averaging %1d FIDs with TE=%1.1f for regression to HSVD basis set\n',numel(index),avgvals(i)); end
        case {'ti','ts'}
            index = find([par.TI] == avgvals(i));
            opt.file_labels{i} = sprintf('TI = %1.0f',avgvals(i)); % overwrite fit result label
            if (opt.verbose), fprintf(1,'Averaging %1d FIDs with %s=%1.1f for regression to HSVD basis set\n',numel(index),opt.avg_group,avgvals(i)); end
        case 'series'
            index = find([par.seriesnum] == avgvals(i));
            opt.file_labels{i} = sprintf('Series = %1d',avgvals(i)); % overwrite fit result label
            if (opt.verbose), fprintf(1,'Averaging %1d FIDs from Series #%1d for regression to HSVD basis set\n',numel(index),avgvals(i)); end
        case {'every','every...'}
            index = (1:opt.fit_avgfids) + (i-1)*opt.fit_avgfids;
            index = index(index <= nfids);
            if (opt.verbose), fprintf(1,'Block averaging FIDs %1d to %1d\n',index(1),index(end)); end
    end
    res(i).par = par(index(1));                     % results will hold params of first FID in averging block
    res(i).par.actual_nex = sum([par(index).nex]);  % if we averaged some FIDs

    % --- Compute FID average ---
    if (opt.verbose)
        fprintf(1,'Fitting FIDs:\n')
        for j=1:numel(index)
            [path,fname,ext] = fileparts(opt.filenames{index(j)});
            [~,tpath] = fileparts(path);
            fprintf(1,'    %s\n',['..' filesep() tpath filesep() fname ext]);
        end
    end
    fid = mean(fids(:,index),2);

    % --- Line broaden  ---
    fidraw = fid;  % need this for noise estimate, un-linebroadened
    if (opt.lb > 0.), fid = fid.*apod; end
    
    % --- calc fit with HSVD basis ---
    if (opt.verbose), fprintf(1,'Regressing FID to %1d component amplitudes\n',comps.ncomps); end
    if (opt.verbose) && (opt.fit_skippts > 0), fprintf(1,'    Ignoring %1d initial FID points in regression\n',opt.fit_skippts); end
    %    amp            = comps.basis((opt.fit_skippts+1):end,:) \ fid((opt.fit_skippts+1):end);
    amp            = comps.basis_inv * fid((opt.fit_skippts+1):end);  % could use backslash, e.g. amp = comps.basis\fid (but not sign. different)
    fitfid         = comps.basis*amp; % This recovers skippts not considered in regression
    comps.amp(:,i) = amp;
    
    % --- compute noise estimate (valid ONLY for DFS spectra!!) ---
    specdiff  = m1dfft(fidraw-fitfid);
    t1        = round(6.5*np/10);   % need points NOT near edge of spectrum where Siemens filtering cuts the noise amp
    t2        = round(8.5*np/10);   % this region is hardcoded for upfield region, which in DFS is mostly signal free
    x         = t1:t2; 
    y         = real(specdiff(x))';
    wstate = warning();
    warning('off','all');               % suppress warnings from polyfit about fit order
        fitorder = 3;
        [p,S] = polyfit(x,y,fitorder);  % de-trend
    warning(wstate);
    trendfit  = polyval(p,x,S);  
    res(i).noise_stdev = std(y - trendfit); % ./ sqrt(numel(x)); 

    % --- Save params of the fitted peaks in results struct ---
    for j=1:peaks.count
        if (~comps.missing_peaks(j))
            p = find(comps.peak_assignment == j);           % these components assigned to peak #j
            res(i).ampcmplx(j) = sum(comps.amp(p,i));       % complex sum of all component amps
            res(i).amp(j)      = abs(res(i).ampcmplx(j));
            if (numel(p) == 1)                              % single component for this peak, so we know lw and freq
                res(i).freq_ppm(j)   = comps.freq_ppm(p(1));
                res(i).lw(j)         = comps.lw(p(1));
                res(i).peakheight(j) = np*2*res(i).amp(j)/pi/res(i).lw(j); % to compute peak height from amp & LW, see https://mathworld.wolfram.com/LorentzianFunction.html
                res(i).snr(j)        = res(i).peakheight(j)/res(i).noise_stdev;  % the "np" term above is the FT scaling factor, since our "amps" are FID(0) values
            else                                            % multiple components for this peak, so we don't know lw and freq
                res(i).freq_ppm(j)   = -99;
                res(i).lw(j)         = -1;
                res(i).peakheight(j) = -1;
                res(i).snr(j)        = -1;
            end
            res(i).ncomps(j) = numel(p);
        else                                                % peak not found, label it w/ amp = -1
            res(i).ampcmplx(j)   = -1;
            res(i).amp(j)        = -1;
            res(i).freq_ppm(j)   = -99;
            res(i).lw(j)         = -1;
            res(i).ncomps(j)     = -1;
            res(i).peakheight(j) = -1;
            res(i).snr(j)        = -1;
        end
    end
    
    % --- Fill in remaining results ---
    res(i).exptype         = opt.exptype;
    res(i).data_name       = opt.data_name;
    res(i).peak_labels     = peaks.labels;
    res(i).status          = 1;  % success

    % --- show a plot of our attempt to measure the spectral noise amp ---
    if (opt.debug && ~isequal(lower(opt.exptype),'water'))
        figure()
        specraw = m1dfft(fidraw);
        specfit = m1dfft(fitfid);
        subplot(2,1,1)
        plot(real(specraw))
        hold on
        plot(real(specfit))
        plot(real(specdiff))       
        hold off
        title('noise estimate')
        axis('tight')

        subplot(2,1,2)
        plot(x,y)
        hold on
        plot(x,trendfit)
        plot(x,y-trendfit)
        hold off
        title(sprintf('sigma = %1.3f',res(i).noise_stdev))
        axis('tight')
    end
    
    % --- Save data for return ---
    if (savedata)
        if (i == 1)
            np            = size(fid,1);
            data.fid      = complex(zeros(np,nout),zeros(np,nout));
            data.fidfit   = complex(zeros(np,nout),zeros(np,nout));
            data.files    = opt.filenames;
            data.peaks    = peaks;
            data.ppm      = ppm;
            data.hz       = hz;
            data.bw       = bw;
            data.f0       = par.f0;
            data.ppmshift = ppmshift;
        end
        data.fid(:,i)    = fid;
        data.fidfit(:,i) = fitfid;
        if (i == nout)
            data.comps    = comps;
            data.status   = 1;
        end
    end

    % --- Get LW of water (need to spline interp to get pseudo-better resolution) ---
    if (isequal(lower(opt.exptype),'water'))
        [res(i).fwhm,s,s_interp,hz_interp,pfwhm1,pfwhm2] = fwhm(m1dfft(fid),hz);
        if (savedata)
            if (i == 1)
                data.water_spec        = zeros(np,nout);
                data.water_spec_interp = zeros(numel(s_interp),nout);
                data.hz_interp         = hz_interp;
                data.pfwhm1(i)         = pfwhm1;
                data.pfwhm2(i)         = pfwhm2;
            end
            data.water_spec(:,i)        = s;
            data.water_spec_interp(:,i) = s_interp;
        end
    end
end
fprintf(1,'-------------------------------------------------------------\n');
fprintf(1,'Done %s()\n',subroutine_name);
fprintf(1,'-------------------------------------------------------------\n\n');
end

% ------------------------------------------------------------------------------------------------
% --- PRINT HSVD COMPONENT PARAMS ---
% ------------------------------------------------------------------------------------------------
function dump_comps(index,freq_ppm,amp,lw,phase,peak_assignments)
ncomps = numel(index);
if (nargin < 6)
    peak_assignments = cell(1,ncomps);
    peak_assignments(:) = {' '};    % no peak assignments provided
end
amp_rel = abs(amp)/max(abs(amp(:)));
disp('comp     freq(ppm)  amp(au)    lw(hz)   phi(deg)       peak');
for k=1:ncomps
    fprintf(1,'%2d) %10.2f %10.4f %10.2f %10.2f %15s\n', index(k),freq_ppm(k),amp_rel(k),lw(k),phase(k),peak_assignments{k});
end
end


% ------------------------------------------------------------------------------------------------
% CHECK OPTIONS STRUCTURE
% ------------------------------------------------------------------------------------------------
function opt = check_options(opt)

opt.fit_order        = checkstruct(opt,'fit_order',[]);             % model order for HSVD
opt.single_component = checkstruct(opt,'single_component','');      % restrict peaks to a single component (can be 'ppm' or 'amp')
opt.missing_peaks_ok = checkstruct(opt,'missing_peaks_ok',0);       % OK if some peaks are not found
opt.subtract_pairs   = checkstruct(opt,'subtract_pairs',0);         % subtract pairs of FIDs before fitting (for spectral editing)
opt.peaktable_shift  = checkstruct(opt,'peaktable_shift',0);        % tweak peak positions in peak table by this global shift
opt.ppm_center       = checkstruct(opt,'ppm_center',4.7);           % make center of ppm scale have this value
opt.center_ppm       = checkstruct(opt,'center_ppm','center');      % find the center of the ppm scale either by 'center' (or 'max')
opt.use_ppm          = checkstruct(opt,'use_ppm',[]);               % use pre-calcualted ppm scale

opt.hsvd_order       = checkstruct(opt,'hsvd_order',60);
opt.model_algo       = checkstruct(opt,'model_algo','HSVD');        % HSVD, HTLS, or HLSVD
opt.use_htls         = checkstruct(opt,'use_htls',0);               % Use HTLS instead of HSVD
opt.separate_model   = checkstruct(opt,'separate_model',0);         % Model blocks of the FIDs separately
opt.sort_by          = checkstruct(opt,'sort_by','none');
opt.hsvd_np          = checkstruct(opt,'hsvd_np',0);                % number of FID points to use in HSVD Hankel matrix
opt.model_avgfids    = checkstruct(opt,'model_avgfids',[]);         % perform HSVD decomp on avg of these fids (e.g. 1:4)
opt.fit_avgfids      = checkstruct(opt,'fit_avgfids',0);            % perform LS fit for amps/phases on average of ALL the FIDs
opt.avg_group        = checkstruct(opt,'avg_group','None');         % Average over some group of spectra ('none', 'all', 'series' ...)
opt.lw_min           = checkstruct(opt,'lw_min',5);                 % weed out any components found with lw < this
opt.lw_max           = checkstruct(opt,'lw_max',3000);              % weed out any components found with lw > this
%opt.amp_min          = checkstruct(opt,'amp_min',0.0005);          % weed out any components found with RELATIVE amp < this
opt.amp_min          = checkstruct(opt,'amp_min',0.0);
opt.ppm_min          = checkstruct(opt,'ppm_min',-20);              % weed out any components found with ppm < this

opt.yoke_peakgroups  = checkstruct(opt,'yoke_peakgroups',0);        % yoke together peak groups according to peaktable yoking params

opt.lb               = checkstruct(opt,'lb',5);                     % LB of fid (i.e. the data to be fit)
opt.baseline_lwmin   = checkstruct(opt,'baseline_lwmin',0);

opt.read_skippts     = checkstruct(opt,'read_skippts',0);           % Ignore initial bad data points when READING Dicom FID data
opt.fit_skippts      = checkstruct(opt,'fit_skippts',0);            % Ignore initial bad data points when FITTING FID data
opt.t2star_extrap    = checkstruct(opt,'t2star_extrap',0);          % extrapoloate the skipped FID points using FWHM (i.e. T2*) estimate
opt.fid_freqalign    = checkstruct(opt,'fid_freqalign',0);          % Frequency align FIDs

opt.use_rawdata        = checkstruct(opt,'use_rawdata',0);           % process raw data files instead of dicom files

opt.verbose            = checkstruct(opt,'verbose',0);
opt.debug              = checkstruct(opt,'debug',0);

if (isempty(opt.model_avgfids)), opt.model_avgfids = [1,1]; end
if (numel(opt.model_avgfids) ~= 2),              fprintf(2,'ERROR: model_avgfids option must either be empty or a 2-element array\n'); return; end
if (opt.model_avgfids(1) > opt.model_avgfids(2)), fprintf(2,'ERROR: model_avgfids element1 > element2!\n'); return; end

end

% ------------------------------------------------------------------------------------------------
% CREATE RESULTS STRUCTURE
% ------------------------------------------------------------------------------------------------
function res = create_res_struct(nout)
res.status          = 0;     % return status failed
res.ampcmplx        = [];
res.amp             = [];
res.exptype         = 'unknown';
res.data_name       = 'unknown';
res.peak_labels     = {'unknown'};
res.peak_assignment = [];
res.freq_ppm        = [];
res.lw              = [];
res.plot_yrangeout  = [];
res.fwhm            = -1;
res.fshifts         = [];
res.par             = [];
res.noise_stdev     = -1;
res.peakheight      = [];
res.snr             = [];
res.combine_source  = 'unknown';
res.fittype         = '';
res = repmat(res,nout,1);
end


