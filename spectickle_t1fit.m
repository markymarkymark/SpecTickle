% ------------------------------------------
% --- Fit for T1 from IR/SatR/ProgSat curve ---
% ------------------------------------------
function res = spectickle_t1fit(res,peaks,fittype,nparams)

st = dbstack; subroutine_name = st.name;
fprintf(1,'\n-------------------------------------------------------------\n');
fprintf(1,'Beginning %s()\n',subroutine_name);
fprintf(1,'-------------------------------------------------------------\n');

tmp = [res.par];
switch lower(fittype)
    case 'invrecovery'
        TI = [tmp.TI];
        fitfunc = @irfit;
        timename = 'TI';
    case 'satrecovery'
        TI = [tmp.TI];
        fitfunc = @srfit;
        timename = 'TS';
    case 'progsat'
        TI = [tmp.TR];
        fitfunc = @srfit;
        timename = 'TR';
end
nt     = numel(TI);
if (numel(unique(TI)) < (nparams+1)), res = []; fprintf(2,'ERROR: %s fit needs at least %1d different %s values\n',fittype,(nparams+1),timename); return; end
amps     = [res(:).amp];
amps     = reshape(amps, peaks.count,nt);   % pull our peak fit amps from result struct array
ampfit   = amps * 0;
ampfixed = amps * 0;
TIfine   = linspace(0,max(TI(:)),nt*10);
finefit  = zeros(peaks.count,nt*10);
for i = 1:peaks.count
    [t1,fit,fitparams,rnorm,fitfine,datafixed] = fitfunc(TI,amps(i,:),TIfine,nparams);
    res(1).t1fit_M0(i)     = fitparams(2);   % first res() structure holds fitted params for all nt res() structures
    res(1).t1fit_k(i)      = fitparams(3);
    res(1).t1fit_t1(i)     = t1;
    res(1).t1fit_rnorm(i)  = rnorm;
    res(1).t1fit_rsq(i)    = corr2(datafixed,fit)^2;  % NOTE: the input data MAY have been sorted and/or negated at points
    res(1).TIfine          = TIfine;
    ampfit(i,:)            = fit;
    finefit(i,:)           = fitfine;
    ampfixed(i,:)          = datafixed;
    fprintf('Found T1=%1.1f k=%1.2f\n',fitparams(1),fitparams(3));
end
% return the fitted amps and TIs
for j = 1:nt
    res(j).TI       = TI(j);
    res(j).ampfit   = ampfit(:,j)';
    res(j).ampfixed = ampfixed(:,j)';
end
res(1).finefit = finefit;  % this is a bit of a kluge
res(1).fittype = fittype;
end


