function outfile = spectickle_printresults(topdir,outroot)

outfile = '';

% --- let user pick folder w/ results ---
if ((nargin < 1) || isempty(topdir))
    if (ispref(mfilename())), lastpath = getpref(mfilename(),'lastpath');
    else, lastpath = pwd();
    end
    topdir = uigetdir(lastpath,'Select the folder with the DFS SVS saved results files');
    if (~ischar(topdir)), return; end           % user hit cancel
    setpref(mfilename(),'lastpath',topdir);     % remember for next time
end

% --- Restore variables from files w/ saved results ---
resfile = [topdir filesep() 'Results.mat'];
optfile = [topdir filesep() 'Options.mat'];
if (isfile(resfile)), res = load(resfile); else, fprintf(2,'ERROR: cannot find file "%s"\n',resfile); return; end
if (isfile(optfile)), opt = load(optfile); else, fprintf(2,'ERROR: cannot find file "%s"\n',optfile); return; end

% --- Figure out what was saved ---
version = checkstruct(res,'version',1.0);
switch(version)
    case 1.0   % Old pre-GUI fits
        % par.par_dfs1 = par.par1;
        % par.par_dfs2 = par.par2;
    case 2.0  
        res = res.res;
        opt = opt.opt;
    case 3.0  
        res = res.res;
        opt = opt.opt;
    otherwise
        fprintf(2,'ERROR: unrecognized version "%1.1f" in the saved data files\n',version);
end

% --- Output WATER fit results? (Only if more than one water result e.g. multiple water files or looping a fit param) ---
if (~isempty(res.res_water) && (numel(res.res_water) > 1))
    if (nargin > 1)
        outfile = sprintf('%s%s%s_water.tsv',topdir,filesep(),outroot);
        fp = fopen(outfile,'w');
    else
        fp = 0;
    end
    
    % --- print out the results in TSV format ---
    printsub(1,opt.opt_water,res.res_water,{},{},fp);
    if (fp ~= 0); fclose(fp); end
end

% --- Print DFS results ---
if (~isempty(res.res_dfs1) && (res.res_dfs1(1).status ~= -1))
    
    % --- handle if no water reference scan ---
    if (~isempty(res.res_water) && (res.res_water(1).status == -1))
        fprintf(2,'WARNING: Dummying water results in %s()\n',mfilename);
        res.res_water.peak_labels   = {'Water'};
        res.res_water.amp           = -1000;
        res.res_water.fwhm          = -1;
        opt.opt_water.fit_skippts   = -1;
        opt.opt_water.hsvd_order    = -1;
        opt.opt_water.t2star_extrap = -1;
    end
    
    % --- Save DFS results to a file? ---
    if (nargin > 1)
        outfile = [topdir filesep() outroot '_svs.tsv'];
        fp = fopen(outfile,'w');
    else
        fp = 0;
    end
    
    % --- print out the results in TSV format ---
    printsub(1, opt.opt_dfs1,res.res_dfs1, opt.opt_water(1),res.res_water(1), fp);
    if (~isempty(res.res_dfs2) && (res.res_dfs2(1).status ~= -1))
        printsub(0, opt.opt_dfs2,res.res_dfs2, opt.opt_water(1),res.res_water(1), fp);
    end
    if (fp ~= 0), fclose(fp); end
end

end

% --------------------------------------------------------------------------------
function printsub(print_header, opt,res, opt_water,res_water, fp)

% --- handle contents ---
nfits  = numel(res);
npeaks = numel(res(1).amp);
nruns  = numel(opt);

% --- print tab delimited results labels ---
if (print_header)
    labels = '';
    % --- Acq params ---
    labels = [labels sprintf('%s\t','DataName')];
    labels = [labels sprintf('%s\t','PatName')];
    labels = [labels sprintf('%s\t','PatAge')];
    labels = [labels sprintf('%s\t','ScanDate')];
    labels = [labels sprintf('%s\t','ExpType')];
    labels = [labels sprintf('%s\t','Anatomy')];
    labels = [labels sprintf('%s\t','SeriesNum')];
    labels = [labels sprintf('%s\t','ImageNum')];
    labels = [labels sprintf('%s\t','SeriesName')];
    labels = [labels sprintf('%s\t','VoxDims(mm)')];
    labels = [labels sprintf('%s\t','VoxSize(cc)')];
    labels = [labels sprintf('%s\t','BW(Hz)')];
    labels = [labels sprintf('%s\t','TxRef(V)')];

    labels = [labels sprintf('%s\t','SeqType')];
    labels = [labels sprintf('%s\t','RFPulse')];
    labels = [labels sprintf('%s\t','RefocDur')];
    labels = [labels sprintf('%s\t','SlicePolarity')];
    labels = [labels sprintf('%s\t','CrushPolarity')];
    labels = [labels sprintf('%s\t','CrushCycle')];
    labels = [labels sprintf('%s\t','SampBeforeEcho')];

    labels = [labels sprintf('%s\t','RFwidth(ppm)')];
    labels = [labels sprintf('%s\t','RFcenter(ppm)')];
    labels = [labels sprintf('%s\t','Crusher(mT)')];
    labels = [labels sprintf('%s\t','TI(ms)')];
    labels = [labels sprintf('%s\t','TR(ms)')];
    labels = [labels sprintf('%s\t','TE(ms)')];
    labels = [labels sprintf('%s\t','Navg')];
    % --- Fit params ---
    labels = [labels sprintf('%s\t','RawData')];
    labels = [labels sprintf('%s\t','RawfreqAlign')];
    labels = [labels sprintf('%s\t','RawPrune')];
    labels = [labels sprintf('%s\t','CoilCombine')];
    labels = [labels sprintf('%s\t','CombineSource')];
    labels = [labels sprintf('%s\t','LB(Hz)')];
    labels = [labels sprintf('%s\t','HSVDorder')];
    labels = [labels sprintf('%s\t','ModelPoints')];
    labels = [labels sprintf('%s\t','SkipFitPts')];
    labels = [labels sprintf('%s\t','T2*Extrap')];
    % --- Water Reference results ---
    if (~isempty(res_water))
        labels = [labels sprintf('HSVDorder(%s)\t',res_water.peak_labels{1})];
        labels = [labels sprintf('SkipFitPts(%s)\t',res_water.peak_labels{1})];
        labels = [labels sprintf('T2*Extrap(%s)\t',res_water.peak_labels{1})];
        labels = [labels sprintf('amp(%s)(a.u.)\t',res_water.peak_labels{1})];
        labels = [labels sprintf('FWHM(%s)(Hz)\t',res_water.peak_labels{1})];
    end
    % --- Peak Fit results ---
    labels = [labels sprintf('%s\t','NoiseStdev')];
    for i=1:npeaks, labels = [labels sprintf('amp(%s)(a.u.)\t',res(1).peak_labels{i})];end
    for i=1:npeaks, labels = [labels sprintf('ppm(%s)\t',res(1).peak_labels{i})];end
    for i=1:npeaks, labels = [labels sprintf('lw(%s)(Hz)\t',res(1).peak_labels{i})];end
    for i=1:npeaks, labels = [labels sprintf('ncomps(%s)\t',res(1).peak_labels{i})];end
    for i=1:npeaks, labels = [labels sprintf('SNR(%s)\t',res(1).peak_labels{i})];end
    
    % --- Post-processing fit results ---
    switch lower(res(1).fittype)
        case 't2decay'
            for i=1:npeaks, labels = [labels sprintf('ampfit(%s)(a.u.)\t',res(1).peak_labels{i})]; end
            for i=1:npeaks, labels = [labels sprintf('T2(%s)(ms)\t',res(1).peak_labels{i})]; end
            for i=1:npeaks, labels = [labels sprintf('M0(%s)(a.u.)\t',res(1).peak_labels{i})]; end
            for i=1:npeaks, labels = [labels sprintf('offset(%s)(a.u.)\t',res(1).peak_labels{i})]; end
            for i=1:npeaks, labels = [labels sprintf('r^2(%s)\t',res(1).peak_labels{i})]; end
        case {'invrecovery','satrecovery','progsat'}
            for i=1:npeaks, labels = [labels sprintf('ampfit(%s)(a.u.)\t',res(1).peak_labels{i})]; end
            for i=1:npeaks, labels = [labels sprintf('T1(%s)(ms)\t',res(1).peak_labels{i})]; end
            for i=1:npeaks, labels = [labels sprintf('M0(%s)(a.u.)\t',res(1).peak_labels{i})]; end
            for i=1:npeaks, labels = [labels sprintf('k(%s)(a.u.)\t',res(1).peak_labels{i})]; end
            for i=1:npeaks, labels = [labels sprintf('r^2(%s)\t',res(1).peak_labels{i})]; end
        otherwise
    end   
    disp(labels);
    if (fp ~= 0), fprintf(fp,'%s\n',labels); end
end

% --- print tab delimited results values ---
for j=1:nfits
    m = min([j nruns]);  % could be multiple fits per call to HSVD routine
    vals = '';
    % --- Acq params ---
    vals = [vals sprintf('%s\t',opt(m).data_name)];
    vals = [vals sprintf('%s\t',res(j).par.patname)];
    vals = [vals sprintf('%1d\t',res(j).par.patage)];
    vals = [vals sprintf('%1d\t',res(j).par.scandate)];
    vals = [vals sprintf('%s\t',opt(m).exptype)];
    vals = [vals sprintf('%s\t',opt(m).tissue_type)];
    vals = [vals sprintf('%1d\t',res(j).par.seriesnum)];
    vals = [vals sprintf('%1d\t',res(j).par.imagenum)];
    vals = [vals sprintf('%s\t',res(j).par.seriesname)];
    vals = [vals sprintf('%1.1f,%1.1f,%1.1f\t',res(j).par.voxdims(1),res(j).par.voxdims(2),res(j).par.voxdims(3))];
    vals = [vals sprintf('%1.1f\t',res(j).par.voxsize)];
    vals = [vals sprintf('%1.1f\t',res(j).par.BW)];
    vals = [vals sprintf('%1.0f\t',res(j).par.txrefvolt)];

    vals = [vals sprintf('%s\t',res(j).par.seqtype)];
    vals = [vals sprintf('%s\t',res(j).par.excitepulse)];
    vals = [vals sprintf('%1.1f\t',res(j).par.refocdur)];
    vals = [vals sprintf('%1d\t',res(j).par.slicepolarity)];
    vals = [vals sprintf('%1d\t',res(j).par.crushpolarity)];
    vals = [vals sprintf('%1d\t',res(j).par.crushcycle)];
    vals = [vals sprintf('%1d\t',res(j).par.presamp)];

    vals = [vals sprintf('%1.1f\t',res(j).par.exciteppm)];
    vals = [vals sprintf('%1.1f\t',res(j).par.excitewidth)];
    vals = [vals sprintf('%1.0f\t',res(j).par.crusher)];
    vals = [vals sprintf('%1.1f\t',res(j).par.TI)];
    vals = [vals sprintf('%1.1f\t',res(j).par.TR)];
    vals = [vals sprintf('%1.1f\t',res(j).par.TE)];
    vals = [vals sprintf('%1d\t',res(j).par.actual_nex)];
    % --- Fit params ---
    vals = [vals sprintf('%1d\t',opt(m).use_rawdata)];
    vals = [vals sprintf('%1d\t',opt(m).raw_freqalign)];
    vals = [vals sprintf('%1d\t',opt(m).raw_channelprune)];
    vals = [vals sprintf('%s\t',opt(m).raw_combine_method)];
    vals = [vals sprintf('%s\t',res(j).combine_source)];
    vals = [vals sprintf('%1d\t',opt(m).lb)];
    vals = [vals sprintf('%1d\t',opt(m).hsvd_order)];
    vals = [vals sprintf('%1d\t',opt(m).hsvd_np)];
    vals = [vals sprintf('%1d\t',opt(m).fit_skippts)];
    vals = [vals sprintf('%1d\t',opt(m).t2star_extrap)];
    % --- Water Reference results ---
    if (~isempty(res_water))
        if (j == 1) % only print water results once per data set
            vals = [vals sprintf('%1d\t',opt_water.hsvd_order)];
            vals = [vals sprintf('%1d\t',opt_water.fit_skippts)];
            vals = [vals sprintf('%1d\t',opt_water.t2star_extrap)];
            vals = [vals sprintf('%1.3f\t',res_water.amp)];
            vals = [vals sprintf('%1.1f\t',res_water.fwhm)];
        else
            vals = [vals sprintf('\t')]; % blanks for these columns
            vals = [vals sprintf('\t')];
            vals = [vals sprintf('\t')];
            vals = [vals sprintf('\t')];
            vals = [vals sprintf('\t')];
        end
    end
    % --- Peak Fit results - 
    vals = [vals sprintf('%g\t',res(j).noise_stdev)]; % Noise
    switch opt(m).exptype  % Amp
        case {'Sel_InvRecovery','nonSel_InvRecovery'}
            for i=1:npeaks, vals = [vals sprintf('%1.4f\t',res(j).ampcorrected(i))]; end
        otherwise
            for i=1:npeaks, vals = [vals sprintf('%1.4f\t',res(j).amp(i))]; end
    end    
    for i=1:npeaks  % PPM
        vals = [vals sprintf('%1.2f\t',res(j).freq_ppm(i))];
    end
    for i=1:npeaks  % LW
        vals = [vals sprintf('%1.1f\t',res(j).lw(i))];
    end
    for i=1:npeaks  % Ncomponents
        vals = [vals sprintf('%1d\t',res(j).ncomps(i))];
    end
	for i=1:npeaks, vals = [vals sprintf('%1.2f\t',res(j).snr(i))]; end

    % --- Post-processing fit results ---
    switch lower(res(1).fittype)
        case 't2decay'
            for i=1:npeaks, vals = [vals sprintf('%1.2f\t',res(j).ampfit(i))]; end
            if (j == 1)
                for i=1:npeaks, vals = [vals sprintf('%1.4f\t',res(j).t2fit_t2(i))]; end
                for i=1:npeaks, vals = [vals sprintf('%1.4f\t',res(j).t2fit_M0(i))]; end
                for i=1:npeaks, vals = [vals sprintf('%1.4f\t',res(j).t2fit_offset(i))]; end
                for i=1:npeaks, vals = [vals sprintf('%1.4f\t',res(j).t2fit_rsq(i))];end
            end
        case {'invrecovery','satrecovery','progsat'}
            for i=1:npeaks, vals = [vals sprintf('%1.2f\t',res(j).ampfit(i))]; end
            if (j == 1)
                for i=1:npeaks, vals = [vals sprintf('%1.1f\t',res(j).t1fit_t1(i))]; end
                for i=1:npeaks, vals = [vals sprintf('%1.2f\t',res(j).t1fit_M0(i))]; end
                for i=1:npeaks, vals = [vals sprintf('%1.4f\t',res(j).t1fit_k(i))];end
                for i=1:npeaks, vals = [vals sprintf('%1.4f\t',res(j).t1fit_rsq(i))];end
            end
        otherwise
    end 
    disp(vals);
    if (fp ~= 0), fprintf(fp,'%s\n',vals); end
end
end







