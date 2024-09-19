% -------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------
function spectickle_make_plots(opt_plot,opt,res,data)

fprintf(1,'\n------------------------------------------------------------\n');
fprintf(1,'Beginning %s()\n',mfilename());
fprintf(1,'-------------------------------------------------------------\n');

opt_plot = check_options(opt_plot);

% --- WATER ---
if (data.data_water.status == 1 && opt_plot.plot_waterfit)
    if (opt_plot.verbose), fprintf(1,'Plotting water reference spectra\n'); end
    tmp1 = opt_plot.plot_xrange; opt_plot.plot_xrange = [];  % always plot water fit over default ppm range
    tmp2 = opt_plot.fid_scale;   opt_plot.fid_scale   = 1;   % don't scale water to anything
    plot_engine(opt_plot, opt.opt_water,res.res_water,data.data_water);
    opt_plot.plot_xrange = tmp1;
    opt_plot.fid_scale   = tmp2;
    fprintf(1,'\n')
end

% --- scale DFS FIDs to the water amp ---
if (data.data_water.status == 1)
    opt_plot.fid_scale = 110E3/res.res_water(1).amp; 
    if (opt_plot.verbose), fprintf(1,'Scaling metabolite plots to water reference\n'); end
end 

% --- SVS/DFS ---
if (opt_plot.verbose), fprintf(1,'Plotting metabolite spectra\n'); end
nruns = numel(data.data_dfs1); % if SEPARATE_HSVD=1, then there are multiple data & results structs
if (nruns == 1)
    plot_engine(opt_plot, opt.opt_dfs1,res.res_dfs1,data.data_dfs1);
else
    for i=1:nruns
        optx = opt.opt_dfs1(i);
        optx.file_labels = {optx.file_labels{i}}; % kluge 
        plot_engine(opt_plot, optx ,res.res_dfs1(i),data.data_dfs1(i));
    end
end

% -------------------------------------------------------------------------------------
% --- Post-processing results - T2 or T1 fits ---
% -------------------------------------------------------------------------------------
res   = res.res_dfs1;
peaks = data.data_dfs1.peaks;
switch lower(res(1).fittype)
    case {'satrecovery','progsat','invrecovery'}
        Time     = [res.TI];
        nT       = numel(Time);
        Timefine = res(1).TIfine; 
        xlabel_string = 'TI/TS/TR (ms)';
        for i = 1:peaks.count, title_string{i} = sprintf('%s: T1 = %1.1f (ms)',peaks.labels{i},res(1).t1fit_t1(i)); end
    case 't2decay'
        Time     = [res.TE];
        nT       = numel(Time);
        Timefine = res(1).TEfine; 
        xlabel_string = 'TE (ms)';
        for i = 1:peaks.count, title_string{i} = sprintf('%s: T2 = %1.1f (ms)',peaks.labels{i},res(1).t2fit_t2(i)); end
    otherwise
        return
end
amp    = reshape([res(:).ampfixed],peaks.count,nT);
ampfit = reshape([res(:).ampfit],peaks.count,nT);
for i = 1:peaks.count
    if (i == 1), f = figure(); set(f,'Position',round(f.Position .* [1 0.25 1.5 1.5])); end
    if (peaks.count > 3)                % if more than 3 peaks, split into to columns of plots
        nrow = ceil(peaks.count/2);
        ncol = 2;
    else
        nrow = peaks.count;
        ncol = 1;
    end
    tmp = 1:(nrow*ncol); plotindex = reshape(tmp,ncol,nrow)';  % complicated way to increment subplots by rows frist
    subplot(nrow,ncol,plotindex(i))
    plot(Time,amp(i,:),'bo');
    hold on
    plot(Time,ampfit(i,:),'rs')
    plot(Timefine,res(1).finefit(i,:),'--r')
    hold off
    title(title_string{i},'Interpreter','none')
    xlabel(xlabel_string)
    %set(gca,'ytick',[]);
    axis tight
end
legend('data','fit','Location','SouthEast')
end



% -------------------------------------------------------------------------------------
% --- plot_engine() --
% -------------------------------------------------------------------------------------
function plot_engine(opt_plot, opt,res,data)

% --- Some plot params ---
thick1 = 0.5;
thick2 = 0.5;
thick3 = 0.5;
datacolor = 'k';
pcolors = [[0, 0.4470, 0.7410]
    [0.8500, 0.3250, 0.0980]
    [0.9290, 0.6940, 0.1250]
    [0.4940, 0.1840, 0.5560]
    [0.4660, 0.6740, 0.1880]
    [0.3010, 0.7450, 0.9330]
    [0.6350, 0.0780, 0.1840]
    [0.8500, 0.3250, 0.0980]
    [0.9290, 0.6940, 0.1250]
    [0.4940, 0.1840, 0.5560]
    [0.4660, 0.6740, 0.1880]
    [0.3010, 0.7450, 0.9330]
];

if (opt_plot.verbose), fprintf(1,'\tSpectral scale factor = %g\n',opt_plot.fid_scale); end

% --- peel off some params ---
peaks = data.peaks;
ppm   = data.ppm;
hz    = data.hz;
comps = data.comps;
[np,nfids] = size(data.fid);
[~,ncomps] = size(comps.basis);

% --- Handle zerofilling ---
if ((opt_plot.plot_zerofill > 1) && (~isequal(lower(opt.exptype),'water')))
    if (opt_plot.verbose), fprintf(1,'\tZero filling results for plotting ONLY by a factor of %1d\n',opt_plot.plot_zerofill); end
    % zfill            = zeros(np*(opt_plot.plot_zerofill-1),1); % this is the 1D zerofill block
    % data.modelfid    = [data.modelfid    ; zfill];
    % data.modelfidfit = [data.modelfidfit ; zfill];
    % zfill            = zeros(np*(opt_plot.plot_zerofill-1),nfids);  % this is the 2D zerofill block
    % data.fid         = [data.fid    ; zfill];
    % data.fidfit      = [data.fidfit ; zfill];
    % zfill            = zeros(np*(opt_plot.plot_zerofill-1),ncomps);  % this is the other 2D zerofill block
    % comps.basis      = [comps.basis ; zfill];
    
    data.modelfid    = zerofill_fid(data.modelfid,-opt_plot.plot_zerofill);
    data.modelfidfit = zerofill_fid(data.modelfidfit,-opt_plot.plot_zerofill);
    data.fid         = zerofill_fid(data.fid,-opt_plot.plot_zerofill);
    data.fidfit      = zerofill_fid(data.fidfit,-opt_plot.plot_zerofill);
    [comps.basis,np] = zerofill_fid(comps.basis,-opt_plot.plot_zerofill);

 %   np         = np * opt_plot.plot_zerofill;
    [~,hz,ppm] = spec_xaxis(np,data.bw,data.f0,0,data.ppmshift);
end

% --- Set x-axis plot range ---
if (~isempty(opt_plot.plot_xrange))
    plot_range = opt_plot.plot_xrange;
else
    plot_range = [min(peaks.ppm) - 1., max(peaks.ppm) + 1.];
end
if (plot_range(2) < plot_range(1)), fprintf(2,'ERROR: Plot range must have Min(ppm) < Max(ppm)\n'); return; end
[~,p0] = min(abs(ppm - plot_range(2)));
[~,p1] = min(abs(ppm - plot_range(1)));
ppm0 = ppm(p0);
ppm1 = ppm(p1);
if (opt_plot.verbose), fprintf(1,'\tSetting plot x-axis range to [%1.1f,%1.1f]\n',ppm0,ppm1); end

% ---  Get components in region of interest ---
pdisp = find((comps.freq_ppm >= ppm1) & (comps.freq_ppm <= ppm0)); % only show components within spec range of interest
ndisp = numel(pdisp);

% -----------------------------------
% --- Show initial HSVD fit ---
% -----------------------------------
if (opt_plot.plot_initialfit)
    f = figure(); set(f,'Position',round(f.Position .* [1 0.25 1.5 1.5]))

    fid     = data.modelfid;
    fitfid  = data.modelfidfit;
    spec    = m1dfft(opt_plot.fid_scale * fid);
    fitspec = m1dfft(opt_plot.fid_scale * fitfid);
    dmax    = max(real(spec(p0:p1)));
    dmin    = min(real(spec(p0:p1)));
    drange  = dmax - dmin;
    roffset = dmin - drange/8;    % shift residue plot to 10% below the data plot

    if (opt_plot.plot_fullspec)
        subplot(2,1,1)
        plot(ppm,real(spec),'color',datacolor,'Linewidth',thick3);
        hold on
        plot([ppm(1),ppm(end)],[0,0],'--k')
        hold off
        xlabel('ppm')
        title(sprintf('%s\nEntire Spectrum',opt.file_labels{1}),'Interpreter','none');
        set(gca,'XDIR','reverse')
        if (~opt_plot.plot_yaxis), set(gca,'ytick',[]); end
        set(gca(), 'XLimSpec', 'Tight');

        subplot(2,1,2)
    end

    plot(ppm(p0:p1),real(spec(p0:p1)),'color',datacolor,'Linewidth',thick3);
    hold on
    plot(ppm(p0:p1),real(fitspec(p0:p1)),'color','r','Linewidth',thick1);
    if (opt_plot.plot_residue), plot(ppm(p0:p1),real(spec(p0:p1)-fitspec(p0:p1))+roffset,'color','k','Linewidth',thick1); end
    plot([ppm(p0),ppm(p1)],[0,0],'--k')
    hold off
    xlabel('ppm')
    title('Initial HSVD Fit');
    if (opt_plot.plot_legends); legend({'spectrum','fit','residue'}); end
    set(gca,'XDIR','reverse')
    if (~opt_plot.plot_yaxis), set(gca,'ytick',[]); end
    if (opt_plot.plot_mastertitle), master_plot_title(opt.data_name,12,0); end
    set(gca(), 'XLimSpec', 'Tight');

    %tightfig();
end

% -----------------------------------
% --- Show spec fits for each FID ---
% -----------------------------------
nplotmax = 30;
if (nfids > nplotmax)
    nfids = nplotmax;
    fprintf(2,'NOTE: number of fitted spectra to plot > %1d. Only plotting first %1d\n',nplotmax,nplotmax);
end

% --- phase component plots to largest fitted peak? ---
% if (opt_plot.plot_autophase)
%     [~,ind]  = max([res(1).amp]);              % find biggest fitted peak
%     maxphase = angle(res(1).ampcmplx(ind));    % get it's phase
%     plot_phase = exp(complex(0,-maxphase));
% else
%     plot_phase = 1;
% end

% --- figure out layout of plots ---
total_plots = opt_plot.plot_totalfit + opt_plot.plot_comps + opt_plot.plot_peaks + opt_plot.plot_peaksum;
switch total_plots
    case 4
        m = 2; n = 2; plotorder = [1,3,2,4];
    case 3
        m = 3; n = 1; plotorder = [1,2,3];
    case 2
        m = 2; n = 1; plotorder = [1,2];
    case 1
        m = 1; n = 1; plotorder = 1;
    otherwise
        fprintf(2,'ERROR: Must select at least one of "opt.plot_totalfit", "opt.plot_comps", "opt.plot_peaks" and "opt.plot_peaksum"\n');
        return
end

% --- Run through each FID and it's fit ---
for i=1:nfids
    fid     = data.fid(:,i);
    fidfit  = data.fidfit(:,i);
    spec    = m1dfft(opt_plot.fid_scale * fid);
    specfit = m1dfft(opt_plot.fid_scale * fidfit);

    % --- Build baseline spectrum ---
    p = find(comps.peak_assignment == (peaks.count+1));   % these are unassigned, i.e. baseline components
    baseline = complex(zeros(np,1),zeros(np,1));
    for j=1:numel(p)
        k = p(j);
        baseline = baseline + m1dfft(opt_plot.fid_scale * comps.amp(k,i) * comps.basis(:,k));
    end

    % --- build each peak fitted spectrum and the sum of all peaks ---
    peakspec = complex(zeros(np,peaks.count),zeros(np,peaks.count));
    peaksum  = complex(zeros(np,1),zeros(np,1));
    for j=1:peaks.count
        if (~comps.missing_peaks(j))
            p = find(comps.peak_assignment == j);   % these components are assigned to peak #j
            for k=1:numel(p)
                r = p(k);
                peakspec(:,j) = peakspec(:,j) +  m1dfft(opt_plot.fid_scale * comps.amp(r,i) * comps.basis(:,r));
            end
        end
        peaksum = peaksum + peakspec(:,j);
    end

    % --- set plot scale, and remember it if fixedscale option ---
    if (~isempty(opt_plot.plot_yrange))
        dmin  = opt_plot.plot_yrange(1);
        dmax  = opt_plot.plot_yrange(2);
        Yshift  = dmax - dmin;
    elseif (i == 1) || (opt_plot.plot_fixedscale == 0)  % autoscale or fix the y-axis range
        dmax    = max(real(spec(p0:p1)));
        dmin    = min(real(spec(p0:p1)));
        Yshift  = dmax - dmin;
        Yshift  = Yshift/2;
    end

    if (opt_plot.plot_autophase)
        [~,phi0] = phase_zero(spec(p0:p1));
        plot_phase = exp(complex(0,phi0*pi/180));
    else
        plot_phase = 1;
    end

    % --- initialize figure and title ---
    title_prefix = sprintf('%s\n',opt.file_labels{i});
    f = figure();
    if (total_plots == 4), set(f,'Position',round(f.Position .* [1 0.4 2.0 1.25]));        % 2x2 plot, so make wide
    else,                  set(f,'Position',round(f.Position .* [1 0.4 1.25 1.8])); end    % 3x1, 2x1, or 1x1, so make tall

    % -------------------------------
    % --- 1) plot Spec and total fit ---
    % -------------------------------
    plotnum = 1;
    if (opt_plot.plot_totalfit)
        subplot(m,n,plotorder(plotnum))
        plot(ppm(p0:p1),real(plot_phase*spec(p0:p1)),'color',datacolor,'Linewidth',thick3);
        hold on
        plot(ppm(p0:p1),real(plot_phase*specfit(p0:p1)),'color','r','Linewidth',thick1);
        if (opt_plot.plot_residue)
            residue = real(plot_phase*spec(p0:p1) - plot_phase*specfit(p0:p1));
            plot(ppm(p0:p1),residue + dmin - Yshift*opt_plot.plot_yspacing,'color','k','Linewidth',thick2);
        end
        hold off
        xlabel('ppm')
        set(gca,'XDIR','reverse')
        set(gca(), 'XLimSpec', 'Tight');

        ylims = get(gca(), 'YLim');   % force default y plot range minimum to go 10% lower
        yr    = ylims(2) - ylims(1);
        ylims(1) = ylims(1) - yr/10;
        ylim(ylims);

        if (i == 1), Yrange{plotnum} = get(gca(), 'YLim'); end
        if     (~isempty(opt_plot.plot_yrange)), ylim(opt_plot.plot_yrange);
        elseif (opt_plot.plot_fixedscale),             ylim(Yrange{plotnum});                    end
        if (~opt_plot.plot_yaxis), set(gca,'ytick',[]); end
        
        if (opt_plot.plot_legends), legend({'spectrum','fit','residue'},'Location','northwest'); end
        title_string = sprintf('%s Total Fit',title_prefix);
        title(title_string,'Interpreter','none');
        plotnum = plotnum+1;
    end

    % ------------------------------------
    % --- 2) plot Spec and component fits ---
    % ------------------------------------
    if (opt_plot.plot_comps)
        subplot(m,n,plotorder(plotnum))
        plot(ppm(p0:p1),real(plot_phase*spec(p0:p1)),'color',datacolor,'Linewidth',thick3);
        hold on
        for j=1:ndisp
            compspec = m1dfft(opt_plot.fid_scale * comps.amp(pdisp(j),i) * comps.basis(:,pdisp(j)));
            if (comps.peak_assignment(pdisp(j)) <= peaks.count)
                colorindex = mod(comps.peak_assignment(pdisp(j)),numel(pcolors))+1;  % makes all components of a peak the same color
                linestyle = '-'; linewidth = thick2;             % a peak component
%                plot(ppm(p0:p1),real(plot_phase*compspec(p0:p1)),'color',pcolors(colorindex,:),'LineWidth',linewidth,'LineStyle',linestyle);
                plot(ppm(p0:p1),real(plot_phase*compspec(p0:p1)) + dmin - (j/ndisp)*Yshift*opt_plot.plot_yspacing,'color',pcolors(colorindex,:),'LineWidth',linewidth,'LineStyle',linestyle);
            elseif ((comps.peak_assignment(pdisp(j)) == (peaks.count+1)) && (~opt_plot.plot_totalbaseline)) % plot the baseline
                color = 'k';
                linestyle = '--'; linewidth = thick1;            % a baseline component
                plot(ppm(p0:p1),real(plot_phase*compspec(p0:p1)),'color',color,'LineWidth',linewidth,'LineStyle',linestyle);
            end
        end
        if (opt_plot.plot_totalbaseline)
            color = 'k';
            linestyle = '--'; linewidth = thick1;            % a baseline component
            plot(ppm(p0:p1),real(plot_phase*baseline(p0:p1)),'color',color,'LineWidth',linewidth,'LineStyle',linestyle);
        end
        if (opt_plot.plot_residue)
            residue = real(plot_phase*spec(p0:p1) - plot_phase*specfit(p0:p1));
            plot(ppm(p0:p1),residue + dmin - 1.5*Yshift*opt_plot.plot_yspacing,'color','k','Linewidth',thick2);
        end

        hold off
        xlabel('ppm')
        set(gca(), 'XDIR',     'reverse')
        set(gca(), 'XLimSpec', 'Tight');

        ylims = get(gca(), 'YLim');   % force default y plot range minimum to go 10% lower
        yr    = ylims(2) - ylims(1);
        ylims(1) = ylims(1) - yr/10;
        ylim(ylims);

        if     (i == 1),                         Yrange{plotnum} = get(gca(), 'YLim'); end   % remember Y plot range for later if multiple FIDs to plot w/ same scale
        if     (~isempty(opt_plot.plot_yrange)), ylim(opt_plot.plot_yrange);
        elseif (opt_plot.plot_fixedscale),       ylim(Yrange{plotnum}); end
        if     (~opt_plot.plot_yaxis),           set(gca,'ytick',[]); end

        if (plotnum == 1), title_string = sprintf('%s Peak Components & Baseline',title_prefix);
        else,              title_string = sprintf('Peak Components & Baseline'); end
        title(title_string,'Interpreter','none');
        plotnum = plotnum+1;
    end

    % -----------------------------------------------
    % --- 3) plot (Spec - baseline) & total of peak fits ---
    % -----------------------------------------------
    if (~isempty(opt_plot.plot_yrange))
        drange  = opt_plot.plot_yrange(2) - opt_plot.plot_yrange(1);
        Yshift   = drange;    % shift residue and components below the data plot
    elseif (i == 1) || (opt_plot.plot_fixedscale == 0)  % autoscale or fix the y-axis range
        dmax    = max(real(spec(p0:p1) - baseline(p0:p1)));
        dmin    = min(real(spec(p0:p1) - baseline(p0:p1)));
        drange  = dmax - dmin;
        Yshift  = drange;    % shift residue and components below the data plot
    end

    if (opt_plot.plot_autophase)
        [~,phi0] = phase_zero(spec(p0:p1) - baseline(p0:p1));
        plot_phase = exp(complex(0,phi0*pi/180));
    else
        plot_phase = 1;
    end

    if (opt_plot.plot_peaksum)
        subplot(m,n,plotorder(plotnum))
        plot(ppm(p0:p1),real(plot_phase*(spec(p0:p1) - baseline(p0:p1))),'color',datacolor,'Linewidth',thick3);
        hold on
%        plot(ppm(p0:p1),real(plot_phase*(peaksum(p0:p1) - baseline(p0:p1))) ,'color','r','Linewidth',thick1);
        plot(ppm(p0:p1),real(plot_phase*(peaksum(p0:p1) - 0*baseline(p0:p1))) ,'color','r','Linewidth',thick1);
        if (opt_plot.plot_residue), plot(ppm(p0:p1),real(plot_phase*(spec(p0:p1) - baseline(p0:p1) - peaksum(p0:p1)) - Yshift*opt_plot.plot_yspacing),'color','k','Linewidth',thick2); end   % residue
        hold off
        xlabel('ppm')
        set(gca(),'XDIR','reverse')
        set(gca(),'XLimSpec', 'Tight');

        if (i == 1), Yrange{plotnum} = get(gca(), 'YLim'); end
        if     (~isempty(opt_plot.plot_yrange)), ylim(opt_plot.plot_yrange);
        elseif (opt_plot.plot_fixedscale),             ylim(Yrange{plotnum});                    end
        if (~opt_plot.plot_yaxis), set(gca,'ytick',[]); end

        location = 'northwest';
        if (opt_plot.plot_legends)
            if (opt_plot.plot_residue), legend({'spectrum','fit','residue'},'Location',location); 
            else,                       legend({'spectrum','fit'},'Location',location);  end
        end

        if (opt_plot.plot_autophase), title_substring = 'Sum of Peaks - Baseline Corrected & Phased';
        else,                    title_substring = 'Sum of Peaks - Baseline Corrected'; end
        if (plotnum == 1),       title_string    = sprintf('%s %s',title_prefix,title_substring);
        else,                    title_string    = title_substring; end
        title(title_string,'Interpreter','none');
        plotnum = plotnum+1;
    end

    % -----------------------------------------------
    % --- 4) plot Spec-baseline & individual peak fits ---
    % -----------------------------------------------
    if (opt_plot.plot_peaks)
        subplot(m,n,plotorder(plotnum))
        if (~isequal(lower(opt.exptype),'water'))
            plot(ppm(p0:p1),real(plot_phase*(spec(p0:p1) - baseline(p0:p1))),'color',datacolor,'Linewidth',thick3);
            hold on
            for j=1:peaks.count
                if (~comps.missing_peaks(j))
                    colorindex = mod(j,numel(pcolors))+1;
%                    plot(ppm(p0:p1),real(plot_phase*peakspec(p0:p1,j)) + j*roffsetX*opt_plot.plot_yspacing/5,'color',pcolors(colorindex,:),'Linewidth',thick2);
                    plot(ppm(p0:p1),real(plot_phase*peakspec(p0:p1,j)) + dmin - j*Yshift*opt_plot.plot_yspacing,'color',pcolors(colorindex,:),'Linewidth',thick2);
                end
            end
            hold off

            xlabel('ppm')
            set(gca(),'XDIR','reverse')
            set(gca(),'XLimSpec', 'Tight');

            if (i == 1), Yrange{plotnum} = get(gca(), 'YLim'); end                 
            % if     (~isempty(opt_plot.plot_yrange)), ylim(opt_plot.plot_yrange);  
            % elseif (opt_plot.plot_fixedscale),             ylim(Yrange{plotnum});                    end
            if (opt_plot.plot_fixedscale),             ylim(Yrange{plotnum});                    end
            if (~opt_plot.plot_yaxis), set(gca,'ytick',[]); end

            location = 'northwest';
            if (opt_plot.plot_legends), legend({'spectrum',peaks.labels{~comps.missing_peaks}},'Location',location); end

            if (opt_plot.plot_autophase), title_substring = 'Individual Peaks - Baseline Corrected & Phased';
            else,                    title_substring = 'Individual Peaks - Baseline Corrected'; end
            if (plotnum == 1),       title_string    = sprintf('%s %s',title_prefix,title_substring);
            else,                    title_string    = title_substring; end
            title(title_string,'Interpreter','none');

            % --- plot water FWHM calc instead ---
        else
            plot(hz,data.water_spec(:,i),'-o')
            hold on
            plot(data.hz_interp,data.water_spec_interp(:,i),'r')
            plot([data.hz_interp(data.pfwhm1(i)),data.hz_interp(data.pfwhm2(i))],[data.water_spec_interp(data.pfwhm1(i),i),data.water_spec_interp(data.pfwhm2(i),i)],'-r*','Linewidth',thick2)
            hold off
            axis([hz(p0),hz(p1),-0.1,1.1])
            xlabel('Hz')
            set(gca,'ytick',[])
            title(sprintf('Water (auto-phased) FWHM = %1.1f',res(i).fwhm))
        end
        res(i).plot_yrangeout = get(gca, 'YLim');
        plotnum = plotnum+1;
    end
    set(gcf,'color','w');
    if (opt_plot.plot_mastertitle), master_plot_title(opt.data_name,12,0); end
    set(gca(), 'XLimSpec', 'Tight');
    %tightfig();
end
end


% ------------------------------------------------------------------------------------------------
% CHECK OPTIONS STRUCTURE
% ------------------------------------------------------------------------------------------------
function opt = check_options(opt)
% --- options for plots, labeling data, etc, ... ---
opt.verbose            = checkstruct(opt,'verbose',1);
opt.debug              = checkstruct(opt,'debug',1);
opt.make_plots         = checkstruct(opt,'make_plots',1);
opt.plot_waterfit      = checkstruct(opt,'plot_waterfit',1);
opt.plot_initialfit    = checkstruct(opt,'plot_initialfit',1);
opt.plot_fullspec      = checkstruct(opt,'plot_fullspec',1);
opt.plot_totalfit      = checkstruct(opt,'plot_totalfit',1);
opt.plot_comps         = checkstruct(opt,'plot_comps',1);
opt.plot_peaks         = checkstruct(opt,'plot_peaks',1);
opt.plot_peaksum       = checkstruct(opt,'plot_peaksum',1);
opt.plot_fixedscale    = checkstruct(opt,'plot_fixedscale',0);
opt.plot_totalbaseline = checkstruct(opt,'plot_totalbaseline',1);       % plots the composite baseline, not the individual components
opt.plot_legends       = checkstruct(opt,'plot_legends',1);
opt.plot_yaxis         = checkstruct(opt,'plot_yaxis',0);
opt.data_name          = checkstruct(opt,'data_name','unknown');        % toplevel data name, e.g. subject_date, ...
opt.exptype            = checkstruct(opt,'exptype','unknown');          % data set info, e.g. "sel_inv", "nonsel_inv", ...
opt.file_labels        = checkstruct(opt,'file_labels',{});             % labels for each data file
opt.plot_residue       = checkstruct(opt,'plot_residue',1);
opt.plot_excitewindow  = checkstruct(opt,'plot_excitewindow',0);
opt.plot_autophase     = checkstruct(opt,'plot_autophase',0);
opt.plot_xrange        = checkstruct(opt,'plot_xrange',[]);
opt.plot_yrange        = checkstruct(opt,'plot_yrange',[]);
opt.plot_mastertitle   = checkstruct(opt,'plot_mastertitle',0);
opt.plot_yspacing      = checkstruct(opt,'plot_yspacing',0.5);
opt.plot_zerofill      = checkstruct(opt,'plot_zerofill',1);
opt.fid_scale          = checkstruct(opt,'fid_scale',1);
end
