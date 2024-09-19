function [fids,rawfids,coeffs,pars,t,hz,ppm,fshift,whitemat,chanphase] = siemens_raw_mrs(rawfiles,coeffs,fshift,whitemat,chanphase,opt)
% SIEMENS_RAW_MRS
% Generate coil combined and signal averaged FIDs from Siemens MRS raw data TWIX files
%       NOTE: If multiple raw files are provided, the FIRST will used to determine coil combination coeffs for ALL the rest.
%
% Created: Mark A. Elliott, PhD
%   melliott@upenn.edu
%   https://www.med.upenn.edu/CAMIPM/mark-elliott.html
fids = []; rawfids = []; pars = [];  t = []; hz = []; ppm = []; 
if (nargin < 1)
    rawfiles = uigetfile_plus(mfilename(),'rawdata',pwd(),{'*.dat';'*.*'},'Select TWIX rawdata file');
    if isempty(rawfiles), return; end
end
if (nargin < 2), coeffs   = []; end
if (nargin < 3), fshift   = []; end
if (nargin < 4), whitemat = []; end
if (nargin < 5), chanphase= []; end
if (nargin < 6), opt      = []; end

% --- Check options ---
opt = check_options(opt);

% --- Run through all raw files ---
if (ischar(rawfiles)), rawfiles = {rawfiles}; end
nfiles = numel(rawfiles);
for i=1:nfiles
    [fid,rawfid,coeffs,t,hz,ppm,par,fshift,whitemat,chanphase] = dfs_rawdata_fid(rawfiles{i},coeffs,fshift,whitemat,chanphase,opt);
    if (isempty(fid) && isempty(rawfid)), fids = []; return; end  % ERROR
    if (i == 1) % Make room for results
        [npts,nmeas] = size(fid);
        fids         = complex(zeros(npts,nmeas,nfiles));
        pars         = repmat(par,1,nfiles);  % now struct(nmeas,nfiles)
    end
    if any(size(fid) ~= [npts,nmeas]), fprintf(2,'ERROR: multiple rawdata files but different data sizes!\n'); fids = []; return; end
    fids(:,:,i) = fid;
    pars(:,i)   = par;
end
end

% ------------------------------------------------------------------------------------------------
% DFS_RAWDATA_FID 
%   Process a TWIX SVS file
% ------------------------------------------------------------------------------------------------
function [fid,rawfid,c,t,hz,ppm,par,fshift,whitemat,chanphase] = dfs_rawdata_fid(rawfile,coeffs,fshift,whitemat,chanphase,opt)
fid = []; rawfid = []; c = []; t = []; hz = []; ppm = [];  
[twix,rawfp,twixnoise] = open_rawdata(rawfile,opt);                 % Open raw data file
par = dfs_svs_acqparams(twix,opt);                                  % Get header info
if (~opt.remove_presamp && opt.remove_adcdelay && (par(1).adcdelay < 0))
    fprintf(2,'ERROR: Dont know how to remove ADCdelay points but keep samples-before-echo.\n');
    return
end
[skippts,pt0,droppts]   = get_time_origin(par(1),opt);              % Determine where t=0 is, and how many raw time points to skip/ignore
[rawfid,npts,~,~,nmeas] = read_rawdata(twix,skippts,droppts,opt);   % Read the raw data
fclose(rawfp);                                                      % stupid mapvbvd() leaves the file open!!
[t,hz,ppm]              = get_timefreq(par(1),npts,opt);
if (opt.rawfids_only), return; end                                  % Just return raw data, no averaging or channel combine
[rawfid,fail]           = average_rawdata(rawfid,par(1),coeffs,t,opt); % Handle signal averaging schemes
if (fail), return; end
if (opt.raw_whiten)
    [rawfid,whitemat] = whiten_rawdata(rawfid,whitemat,twixnoise,opt); % Whiten the data
end
if (opt.raw_dephase || ~isempty(chanphase))                         % Do eddy current correction style channel dephasing
    [rawfid,chanphase] = channel_dephase(rawfid,chanphase,opt);
end
if (isempty(coeffs))                                                % Compute coil combine coeffs (and possibly channel frequency shifts)
    [c,fshift] = get_coilcoeffs(rawfid,t,pt0,ppm,opt);
else
    if (opt.verbose), fprintf(1,'    Coil channel combination will use provided coeffs\n'); end
    c = coeffs;                                                     % NOTE that these coeffs are [Nchan x 1], not [Nchan x nfiles]
end
if (~isempty(fshift))                                               % Apply freq shifts to coil channel signals
    rawfid = channel_freqshift(rawfid,fshift,t,opt);
end
for i = 1:nmeas                                                     % combine coil elements
    if (i == 1), fid = complex(zeros(npts,nmeas)); end
    fid(:,i) = squeeze(rawfid(:,:,i)) * c;
end
fid = fid * exp(complex(0,-angle(fid(pt0,1))));                     % Zero constant phase (this is redundant for some methods)
fid = opt.raw_scale * fid;
end

% ------------------------------------------------------------------------------------------------
% NORMALIZE COMPLEX FIDs - either by 1st point mag or by complex amp (i.e. including phase)
% ------------------------------------------------------------------------------------------------
function [fids,amp,phi] = apnorm(rawfids,how)
if (nargin < 2), how = 'amp_phase'; end
[npts,nchan] = size(rawfids);
fids         = complex(zeros(npts,nchan));
switch lower(how)
    case 'amp_phase'
        anorm = rawfids(1,:);
    case 'amp'
        anorm = abs(rawfids(1,:));
end
for k=1:nchan, fids(:,k) = rawfids(:,k) ./ anorm(k); end
amp = abs(anorm).';
phi = angle(anorm).';
end

% ------------------------------------------------------------------------------------------------
% GET START AND END POINTS FOR A RANGE OF PPM 
% ------------------------------------------------------------------------------------------------
function [p0,p1,ppm0,ppm1] = get_ppm_bounds(ppm,ppmcenter,ppmwidth)
ppm0   = ppmcenter - ppmwidth/2;
ppm1   = ppmcenter + ppmwidth/2;
[~,p0] = min(abs(ppm - ppm1));
[~,p1] = min(abs(ppm - ppm0));
end

% ------------------------------------------------------------------------------------------------
% COMPUTE SIMILARITY SCORE FROM RAW FIDS 
% ------------------------------------------------------------------------------------------------
function [d,SIM,order,SIMsorted,m1,m2,m3,p0,p1,data] = compute_similarity(rawfids,ppm,ppmcenter,ppmwidth,opt)
if (nargin < 2), ppm       = []; end
if (nargin < 3), ppmcenter = []; end
if (nargin < 4), ppmwidth  = []; end
if (isempty(ppmcenter)), ppmcenter = 0.0; end
if (isempty(ppmwidth)),  ppmwidth  = 0.4; end

[np,nchan] = size(rawfids);
if (~isempty(ppm))                      % compute only on spectral sub-region
    [p0,p1,ppm0,ppm1] = get_ppm_bounds(ppm,ppmcenter,ppmwidth);
    data              = m1dfft(rawfids);
    data              = data(p0:p1,:);
    if (opt.verbose), fprintf(1,'    Computing Similarity matrix on %1dx%1d spectral data (from %1.1f to %1.1f ppm)\n',p1-p0+1,nchan,ppm0,ppm1); end
else                                    % compute over whole data set
    data = rawfids;
    if (opt.verbose), fprintf(1,'    Computing Similarity matrix on %1dx%1d raw data\n',np,nchan); end
end
SIM = zeros(nchan,nchan);               % compute Similarity matrix
for k=1:nchan
    veck = [real(data(:,k)); imag(data(:,k))];
    for l=k:nchan
        vecl = [real(data(:,l)); imag(data(:,l))];
        SIM(k,l) = corr2(veck,vecl)^2;  % Need to think about R^2 when R<0 ???
        SIM(l,k) = SIM(k,l);
    end
end
d = median(SIM);
[~,order] = sort(d,'descend');
SIMsorted = SIM(order,order);

% --- Some metrics from the Similarity scores ---
x    = 1:nchan;
P    = polyfit(x,d(order),1);
yfit = P(1)*x + P(2);
m1   = -10 * P(1);     % -slope x 10
m2   = (d(order(1)) - d(order(end))) / d(order(1));
yint = P(1)*nchan + P(2);   % the intercept of the line fit to the last channel
pcut = find(d < yint);      % cut channels w/ similarity below the y-intercept
ncut = numel(pcut);
pcut = [order(end-ncut)  pcut]; % add one more low similarity channel to cut
ncut = ncut + 1;
m3   = pcut;

% --- plots ---
if (opt.debug)
    figure()
    subplot(2,1,1)
    imagesc(1-SIMsorted,[0 1])
    xticks(1:numel(d))
    xticklabels(string(order))
    set(gca,'ytick',[]);
    colormap jet
    colorbar
    set(gcf,'color','w');

    subplot(2,1,2)
    plot(x,d(order),'-o')
    hold on;
    plot(x,yfit,'r-.');
        if (opt.plot_prune), plot(x(end-ncut+1:end),d(order(end-ncut+1:end)),'rx'); end
    hold off
    xticks(1:numel(d))
    xticklabels(string(order))
    axis tight
    ylim([0 1]);
    xlabel('channel')
    set(gcf,'color','w');
end
end

% ------------------------------------------------------------------------------------------------
% FREQUENCY SHIFT (AND AMP/PHASE) FITTING OF RAW FIDS USING SIMILARITYSCORE
% ------------------------------------------------------------------------------------------------
function [fshift,amp,phi,m1,m2,m3,fid1,d1] = coilchan_freqalign(rawfids,t,ppm,ppmcenter,ppmwidth,opt)
if (nargin < 3), ppm       = []; end
if (nargin < 4), ppmcenter = []; end
if (nargin < 5), ppmwidth  = []; end

opt.plot_prune = 0;                                                 %  for debug mode plots
[fid0,amp0,phi0] = apnorm(rawfids,'amp_phase');
[~,~,order0] = compute_similarity(fid0,ppm,ppmcenter,ppmwidth,opt); % Compute similarity metric for each coil channel 

[npts,nchan] = size(rawfids);                                       % freq align of each chan to best channel
fid1 = complex(zeros(npts,nchan));   
Nfit    = npts/2;                                                   % number of FID points to use in freq shift fit
skippts = 0;
if (opt.verbose)
    fprintf(1,'    Optimizing frequency shift for %1d coil channels to reference channel %1d\n',nchan-1,order0(1));
    fprintf(1,'       Fitting %1d FID points and skipping %1d initial points\n',Nfit,skippts);
end
fshift = zeros(nchan,1);
amp1   = zeros(nchan,1);
phi1   = zeros(nchan,1);
fidref = fid0(:,order0(1));                                         % use best channel as reference FID
for k=2:nchan
    [fshift(order0(k)),fidfit,pars] = fidfreqshift(t,fid0(:,order0(k)),fidref,[0,1,0],Nfit,skippts);
    amp1(order0(k))   = pars(2);
    phi1(order0(k))   = pars(3);
    fid1(:,order0(k)) = fidfit;
end
amp1(order0(1))   = 1;                                              % these are for the reference FID
phi1(order0(1))   = 0;
fid1(:,order0(1)) = fidref;
if (opt.verbose), fprintf(1,'       Found freq shifts from %1.1f to %1.1f (Hz)\n',min(fshift(:)),max(fshift(:))); end

opt.plot_prune = 1; 
[d1,~,~,~,m1,m2,m3] = compute_similarity(fid1,ppm,ppmcenter,ppmwidth,opt); % Compute Similarity again on freq aligned channels 
amp = (amp0 ./ amp1).';                                             % Return net amp and phase for coil coeffs ---
phi = (phi0 + phi1).';

if (opt.debug), debug_plotwater_beforeafter(fid0,fid1,'Frequency Aligning',opt); end  %  Show before and after plots
end

% ------------------------------------------------------------------------------------------------
% PRUNE BAD COIL CHANNELS BASED ON SIMILARITY METRICS
% ------------------------------------------------------------------------------------------------
function [cout,prunes] = coilchan_prune(c,m1,m2,m3,fids,opt)
global PPM PPM0 PPM1 P0 P1

% --- Remove channels found by "metric 3" ---
cout         = c;
prunes       = m3;
ncut         = numel(prunes);
cout(prunes) = 0;
if (opt.verbose), fprintf(1,'    Pruning %1d coil channels from coil coefficients\n',ncut); end

% --- plot the affects ---
if (opt.debug)
    spec = m1dfft(fids);
    nchan = size(fids,2);

    % --- Show coil channels to be cut v. coil coeff mag ---
    [cmag,corder]   = sort(abs(c).^2,'descend');
    cutindex        = zeros(nchan,1);
    cutindex(m3)    = 1;
    cutindex        = cutindex(corder);
    pcut            = (cutindex == 1);
    x               = 1:nchan;
    figure()
    subplot(2,1,1)
    plot(cmag,'o')
    hold on
    plot(x(pcut),cmag(pcut),'x','Linewidth',2)
    xticks(1:nchan)
    xticklabels(string(corder))
    xlabel('channel')
    title('Coil Combine Coeffs')
    axis tight

    subplot(2,1,2)
    plot(PPM(P0:P1),real(spec(P0:P1,corder(1))),'Linewidth',1)
    hold on
    for k=2:nchan
        if (~any(m3 == corder(k)))
            plot(PPM(P0:P1),real(spec(P0:P1,corder(k))));
        end
    end
    hold off
    xlabel('ppm')
    axis tight
    set(gca,'XDIR','reverse')
    set(gca,'ytick',[],'Ycolor','w','box','off')
    set(gca,'XAxisLocation','bottom', 'box','off')
    set(gcf,'color','w');
    title('After Pruning')
end
end

% ------------------------------------------------------------------------------------------------
% OPEN RAWDATA FILE
% ------------------------------------------------------------------------------------------------
function [twix,rawfp,twixnoise] = open_rawdata(rawfile,opt)
if (opt.verbose), fprintf(1,'Reading raw TWIX file: %s\n',rawfile); end
[twix,rawfp] = me_mapVBVD(rawfile,'verbose',false); % Read raw data
twixnoise = [];                                     % Look for noise data in TWIX
if (iscell(twix))                 
    if (isfield(twix{1},'noise'))
        twixnoise = twix{1}.noise;
    end
    twix = twix{2};                                 % now point to the actual data
elseif (isfield(twix,'noise'))
    twixnoise = twix.noise;
end
end

% ------------------------------------------------------------------------------------------------
% READ RAW DATA
% ------------------------------------------------------------------------------------------------
function [data,npts,nchan,navg,nmeas] = read_rawdata(twix,skippts,droppts,opt)
do_average                 = isequal(lower(opt.raw_avg_method),'simple');   % simple block averging of all NEX
twix.image.flagDoAverage   = do_average;                                    % fastest/easiest to do this on the read in mpaVBVD()
twix.image.flagAverageSets = do_average;                                    % this is for the EJA sequences (they put averages in the 'set' location!!)
twix.image.flagRemoveOS    = opt.remove_oversamp;
data                       = twix.image();
data                       = conj(data);                                    % imag channel needs negation (needed to reverse spectra properly)
npts     = size(data,1);
nchan    = twix.image.NCha;
nmeas    = twix.image.NEco;
navg     = twix.image.NAve;                                  
navg_eja = twix.image.NSet;                                                 % EJA seq uses 'Set' for averages
if (opt.verbose && nmeas > 1), fprintf(2,'\tNOTE: Found multiple measurements (%1d) in rawdata file\n',nmeas); end
if (opt.verbose), fprintf(1,'\tInitial data size [npts,nchan,navg,nmeas] = [%1d,%1d,%1d,%1d]\n',npts,nchan,navg*navg_eja,nmeas); end
if (do_average)  % averages were removed on read()
    if (opt.verbose), fprintf(1,'\tPerforming averaging on rawdata read\n'); end
    navg = 1; navg_eja = 1; 
end                                
data = squeeze(data);
if (navg_eja > 1)                                                           % need to move average before measurements
    data = reshape(data,npts,nchan,nmeas,navg_eja);
    data = permute(data,[1,2,4,3]);
else
    data = reshape(data,npts,nchan,navg,nmeas);                             % averages in right place, just leave singleton
end
if (skippts > 0 || droppts > 0)                                             % discard leading/trailing samples to ignore
    data = data((1+skippts):(npts-droppts),:,:,:);  
end
[npts,nchan,navg,nmeas] = size(data);                                       % returns final dims
if (opt.verbose), fprintf(1,'\tFinal data size [npts,nchan,navg,nmeas] = [%1d,%1d,%1d,%1d]\n',npts,nchan,navg*navg_eja,nmeas); end
end

% ------------------------------------------------------------------------------------------------
% CALCULATE TIME, HZ, and PPM SCALES
% ------------------------------------------------------------------------------------------------
function [t,hz,ppm]  = get_timefreq(par,npts,opt)
global PPM PPM0 PPM1 P0 P1

[t,hz,ppm] = spec_xaxis(npts,par.BW,par.f0,par.adcdelay/1E6,0);

if (opt.debug)          % hang on to time/freq axes for debug routines
    ppmcenter = 0.0;
    ppmwidth  = 0.6;
    PPM       = ppm + 0.0;
    [P0,P1,PPM0,PPM1] = get_ppm_bounds(PPM,ppmcenter,ppmwidth);
end
end

% ------------------------------------------------------------------------------------------------
% AVERAGE RAW UNCOMBINED DATA
% ------------------------------------------------------------------------------------------------
function [rawfids,fail] = average_rawdata(rawfids,par,coeffs,t,opt)
fail = 0;
if (opt.verbose), fprintf(1,'    Averaging FIDs using the "%s" method\n',opt.raw_avg_method); end
[npts,nchan,navg] = size(rawfids);
switch lower(opt.raw_avg_method)
    case 'simple'  % this was accomplished in the Twix read
        % rawfids = mean(rawfids,3); % if not done on TWIX read
    case 'svd'
        % [~,ind] = max(abs(coeffs)); % find channel with largest signal weighting
        % coildat = squeeze(rawfids(:,ind,:));
        % coildat = reshape(coildat,npts,par.phcyclesteps,navg/par.phcyclesteps);  % Average over N-step phase cycle
        % coildat = squeeze(mean(coildat,2));
        % [U,S,V] = svd(coildat,'econ');
        % c       = V(:,1) / S(1,1);  % NOT FINISHED !!
        fprintf(2,'ERROR: SVD averaging scheme is not completely coded yet\n');
        fail = 1;
    case 'block'
        blocksize = opt.raw_avg_block;  
        if (mod(navg,blocksize) ~= 0), fprintf(2,'ERROR: Number of rawdata averages (%1d) is not divisible by blocksize (%1d)\n',navg,blocksize); fail= 1; return; end
        nblocks = navg/blocksize;
        rawfids = reshape(rawfids,npts,nchan,blocksize,nblocks);
        rawfids = squeeze(mean(rawfids,3));  % Average over blocks
        rawfids = rawfids(:,:,1); % KLUGE for now - return only first block
    case 'freq align'
        if (isempty(coeffs)), fprintf(2,'ERROR: %s averaging method requires coil combine coefficients to be provided\n',opt.raw_avg_method); fail = 1; return; end
        blocksize = max([par.phcyclesteps 8]);   % phase cycling block size, but not less than 8
        nblocks = navg/blocksize;
        rawfids = reshape(rawfids,npts,nchan,blocksize,nblocks);
        rawfids = squeeze(mean(rawfids,3));  % Average over N-step phase cycle blocks
        rawfidsX = complex(zeros(npts,nblocks));
        for j=1:nblocks
            rawfidsX(:,j) = squeeze(rawfids(:,:,j)) * coeffs; % combine channels but keep signal averaging blocks
        end
        scale = abs(rawfidsX(1,1));
        rawfidsX = rawfidsX / scale;  % curvefit doesn't do well on really small numbers
        fidref = rawfidsX(:,1);
        fshift = zeros(nblocks,1);
        for j=2:nblocks
            [fshift(j),ffit,fpars] = fidfreqshift(t,rawfidsX(:,j),fidref,[0,1,0],npts/2);
            %if (opt.verbose), fprintf(1,'Frequency drift correction of %5.1f Hz applied to averaging block %2d\n',fshift(j),j); end
        end
        if (opt.verbose), fprintf(1,'    Frequency drift correction applied to %1d signal-average blocks (mean/max/min = %1.1f/%1.1f/%1.1f)\n',nblocks,mean(fshift),min(fshift),max(fshift)); end
        for k=2:nblocks             % now apply to fshifts to uncombined avergaing blocks
            rawfids(:,:,k) = exp(-1i*fshift(k)*2*pi*t) .* rawfids(:,:,k);
        end
        rawfids  = mean(rawfids,3);  % finish averaging - these are coil uncombined!!
        if (opt.debug)
            figure()
            plot(fshift)
            xlabel('averging block')
            ylabel('freq shift (Hz)')
            title(sprintf('Frequency Drift "%s"',par.seriesname),'Interpreter','none')
            axis tight
        end
    otherwise
        fprintf(2,'ERROR: Unknown averaging scheme "%s"\n',opt.raw_avg_method);
        fail = 1;
end
end

% ------------------------------------------------------------------------------------------------
% COMPUTE WHITENING MATRIX
% ------------------------------------------------------------------------------------------------
function whitemat = compute_whitemat(twixnoise,opt)

if (opt.verbose), fprintf(1,'    Computing whitening matrix from coil noise data\n'); end
twixnoise.flagDoAverage = 1;
twixnoise.flagRemoveOS  = 0;
noise                   = twixnoise();    % reads noise w/ signal averaging done, but channels separate
[~,~,nlin]              = size(noise);
if (nlin > 1), noise = mean(noise,3); end % nLIN are really averages?
noise                   = conj(noise);    % imaginary channel needs negation?
npnoise                 = size(noise,1);
% mu_X   = mean(noise, 1);               % (use Laurens van der Maaten, Delft University of Technology)
% X      = bsxfun(@minus, noise, mu_X);  % this is the de-meaned noise data
% Xwhite = X / sqrtm(cov(X));            % this is whitened noise data
% W1     = X \ Xwhite;                   % whitening coeffs
% W1     = W1 / norm(W1);                % make whitening matrix unit scaling
psi = (1/(npnoise-1))*(noise' * noise);  % (use Hansen http://hansenms.github.io/sunrise/sunrise2014/Hansen_ImageReconstruction_ParallelImaging_20140428.pdf)
L   = chol(psi,'lower');
W2  = inv(L);                            % Whitening coeffs
W2  = W2 / norm(W2);                     % make whitening matrix have unit scaling
whitemat = W2;
if (opt.debug)
    figure()
    % imagesc(abs(whitemat))
    % colormap jet
    % title('Whitening Matrix');
    imagesc(abs(psi))
    colormap hot
    title('Coil Noise Correlation Matrix');
end
end

% ------------------------------------------------------------------------------------------------
% WHITEN RAW DATA
% ------------------------------------------------------------------------------------------------
function [rawfids,whitemat] = whiten_rawdata(rawfids,whitemat,twixnoise,opt)
if (~isempty(whitemat))
    if (opt.verbose), fprintf(1,'    Whitening raw data w/ provided whitening matrix\n'); end
    rawfids = (whitemat * rawfids.').';
elseif (isempty(whitemat) && ~isempty(twixnoise))
    whitemat = compute_whitemat(twixnoise,opt);
    if (opt.verbose), fprintf(1,'    Whitening raw data from coil noise data\n'); end
    rawfids = (whitemat * rawfids.').';
else
    if (opt.verbose), fprintf(2,'    WARNING: whitening option ON but no whitening matrix provided or coil noise found\n'); end
end
end

% ------------------------------------------------------------------------------------------------
% APPLY FREQ SHIFTS TO COIL CHANNEL SIGNALS
% ------------------------------------------------------------------------------------------------
function rawfid = channel_freqshift(rawfid,fshift,t,opt)
if (opt.verbose), fprintf(1,'    Applying frequency shift corrections to coil channel data\n'); end
nchan = size(rawfid,2);
for k=1:nchan, rawfid(:,k) = exp(-1i*fshift(k)*2*pi*t) .* rawfid(:,k); end
end

% ------------------------------------------------------------------------------------------------
% GET (water!) COIL CHANNEL DE-PHASING VECTORS
% ------------------------------------------------------------------------------------------------
function chanphase = get_channel_dephasers(rawfid,opt)
if (opt.verbose), fprintf(1,'    Calculating coil channel phase vectors from raw data\n'); end
chanphase = exp(-1i*angle(rawfid));
if (opt.debug), debug_plotwater_beforeafter(rawfid,rawfid .* chanphase,'Water de-Phasing',opt); end
end

% ------------------------------------------------------------------------------------------------
% APPLY COIL CHANNEL DE-PHASING
% ------------------------------------------------------------------------------------------------
function [rawfid,chanphase] = channel_dephase(rawfid,chanphase,opt)
    if (isempty(chanphase)), chanphase = get_channel_dephasers(rawfid,opt); end
    if (opt.verbose), fprintf(1,'    Applying coil channel de-phasing\n'); end
    rawfid = rawfid .* chanphase;
end

% ------------------------------------------------------------------------------------------------
% CALCULATE COIL COMBINE COEFFICIENTS
% ------------------------------------------------------------------------------------------------
function [c,fshift] = get_coilcoeffs(rawfids,t,pt0,ppm,opt)

fshift = [];
[~,nchan,nmeas] = size(rawfids);
if (nmeas > 1)
    fprintf(1,'    Multiple (%1d) measurements in rawdata for coil combining. Using 1st measurement\n',nmeas); 
    rawfids = rawfids(:,:,1);
end
% --- Determine Coil coeffs with fit for freq shift, 0-phase, and amplitude ---
if (opt.raw_freqalign)      
    [fshift,amp1,phi1,m1,m2,m3,fidx] = coilchan_freqalign(rawfids(pt0:end,:),t,ppm,[],[],opt);
    if (opt.raw_ignore_fshifts)
        if (opt.verbose), fprintf(1,'    Discarding freq shifts\n'); end
        fshift = [];     % Experimental - show effect of just amp/phase fitting
    end
end
% --- Use amps and phases from "conventional" FID(0) point ---
if (~opt.raw_freqalign || opt.raw_fshiftonly)
    if (opt.verbose), fprintf(1,'    Determining channel amps/phases from FID(t=0) point\n'); end
    amp1 = abs(rawfids(pt0,:));
    phi1 = angle(rawfids(pt0,:));
else
    if (opt.verbose), fprintf(1,'    Using channel amps/phases from freq align fitting\n'); end
end
% --- Coil coeffs weighting scheme ---
if (opt.verbose), fprintf(1,'    Determining coil combine coeffs w/ "%s" algorithm\n',opt.raw_combine_method); end
switch lower(opt.raw_combine_method)
    case 'signal'
        c = amp1 .* exp(-1i*phi1);                      % |fid(1)| w/ phase negated
        c = c.';
    case {'signalsquared','signal^2'}
        c = (amp1.^2) .* exp(complex(0,-phi1));         % |fid(1)|^2 w/ phase negated
        c = c.';
    case 'svd'                                          % from WSVD (Magnetic Resonance in Medicine 63:881–891 (2010))
        [U,S,V] = svd(rawfids,'econ');
        c       = V(:,1) / S(1,1);                      % Use first principal component as best estimate of coil combines
        % i.e. rawfids * c = U(:,1)
    case {'maxchan','max channel'}
        [~,ind] = max(amp1(:));                         % Use only the coil channel with max signal
        c = zeros(nchan,1);
        c(ind) = 1.0 .* exp(complex(0,-phi1(ind)));
    otherwise
        fprintf(2,'ERROR: Unknown coil combine mode "%s"\n',opt.raw_combine_method);
        fids = [];
        return
end
% --- Prune bad channels from freq align step ---
if (opt.raw_freqalign && opt.raw_channelprune)
    c = coilchan_prune(c,m1,m2,m3,fidx,opt);
end
% --- % normalize so sum(|c|^2) = 1  ---
sumsq = sum(c .* conj(c));
c     = c / sqrt(sumsq);
end

% ------------------------------------------------------------------------------------------------
% DETERMINE DATA TIME ORIGIN
% ------------------------------------------------------------------------------------------------
function [skippts,pt0,droppts] = get_time_origin(par,opt)
skippts      = 0;
adcdelay_pts = 0;
if (par.adcdelay < 0),     adcdelay_pts = -round(par.adcdelay/1e6*par.BW); end % number of FID points in ADCdelay
if (opt.remove_presamp),   skippts = par.presamp; end                          % throw out samples before echo
if (opt.remove_adcdelay),  skippts = skippts + adcdelay_pts; end               % throw out early ADC points
pt0     = par.presamp + adcdelay_pts - skippts + 1;                            % this is the index of t=0 time point
droppts = par.postsamp;                                                        % extra points at the end

if (opt.verbose)
    if (opt.remove_oversamp),  fprintf(1,'    Removing Over Sampling points\n'); end
    if (par.postsamp > 0) || 1,fprintf(1,'    Found %1d samples-after-fid\n',droppts); end
    if (par.presamp > 0) || 1, fprintf(1,'    Found %1d samples-before-echo\n',par.presamp); end
                               fprintf(1,'    Found ADC delay of %1.3f msec (%1d samples)\n',par.adcdelay/1000,adcdelay_pts);
    if (skippts > 0) || 1,     fprintf(1,'    Omitting %1d initial raw FID points\n',skippts); end
                               fprintf(1,'    Using data point %1d for FID t=0\n',pt0);
end
end

% ------------------------------------------------------------------------------------------------
% CHECK OPTIONS STRUCTURE
% ------------------------------------------------------------------------------------------------
function opt = check_options(opt)
if (isempty(opt)), opt.dummy = 1; end
opt.verbose            = checkstruct(opt,'verbose',1);
opt.debug              = checkstruct(opt,'debug',0);
opt.raw_scale          = checkstruct(opt,'raw_scale',1);                 % scale raw data by this factor
opt.rawfids_only       = checkstruct(opt,'rawfids_only',0);              % return only raw data, no averaging or channel combine
opt.raw_combine_method = checkstruct(opt,'raw_combine_method','signal'); % method for determine coil combine coeffs
opt.remove_presamp     = checkstruct(opt,'remove_presamp',1);            % skip initial FID points due to samples-before-echo
opt.remove_oversamp    = checkstruct(opt,'remove_oversamp',0);           % remove OS
opt.remove_adcdelay    = checkstruct(opt,'remove_adcdelay',1);           % Remove ADCdelay points (if negative)
opt.raw_whiten         = checkstruct(opt,'raw_whiten',0);                % whiten coil element signals 
opt.raw_avg_method     = checkstruct(opt,'raw_avg_method','simple');     % averaging mode 
opt.raw_avg_block      = checkstruct(opt,'raw_avg_block',32);            % averaging block size (if 'raw_avg_method' = 'block')
opt.raw_freqalign      = checkstruct(opt,'raw_freqalign',0);             % correct for coil element freq i.e. B0 differences
opt.raw_channelprune   = checkstruct(opt,'raw_channelprune',0);          % delete coil elements w/ significant non-similarity
opt.raw_fshiftonly     = checkstruct(opt,'raw_fshiftonly',0);            % only use fshifts, not additional amp/phase from freq-align fit
opt.raw_ignore_fshifts = checkstruct(opt,'raw_ignore_fshifts',0);        % only use amp/phase from freq-align fit (not fshifts)
opt.raw_dephase        = checkstruct(opt,'raw_dephase',0);               % perform eddy current style dephasing (only for WATER)
end


% ------------------------------------------------------------------------------------------------
% DEBUGGING ROUTINES
% ------------------------------------------------------------------------------------------------

% ------------------------------------------------------------------------------------------------
% ------------------------------------------------------------------------------------------------
function debug_plotwater_beforeafter(rawfids1,rawfids2,titstring,opt)
global PPM PPM0 PPM1 P0 P1
rawfids1 = apnorm(rawfids1);
rawfids2 = apnorm(rawfids2);
spec1    = m1dfft(rawfids1);
spec2    = m1dfft(rawfids2);
nchan    = size(spec2,2);

figure()
subplot(2,1,1)
    plot(PPM(P0:P1),real(spec1(P0:P1,1)))
    hold on
    for k=2:nchan, plot(PPM(P0:P1),real(spec1(P0:P1,k))); end
    hold off
    debug_specplot_options();
    title([titstring ' - BEFORE'])
subplot(2,1,2)
    plot(PPM(P0:P1),real(spec2(P0:P1,1)))
    hold on
    for k=2:nchan, plot(PPM(P0:P1),real(spec2(P0:P1,k))); end
    hold off
    debug_specplot_options();
    title([titstring ' - AFTER'])
end

% ------------------------------------------------------------------------------------------------
% ------------------------------------------------------------------------------------------------
function debug_specplot_options()
xlabel('ppm')
axis tight
hold off
set(gca,'XDIR','reverse')
set(gca,'Xticklabel',[])
set(gca,'ytick',[],'Ycolor','w','box','off')
set(gca,'XAxisLocation','bottom', 'box','off')
set(gcf,'color','w');
end