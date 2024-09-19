% ------------------------------------------
% --- Fit for T1 from inversion or ssaturation recovery curve ---
% ------------------------------------------
function res = spectickle_t2fit(res,peaks,nparams)

st = dbstack; subroutine_name = st.name;
fprintf(1,'\n-------------------------------------------------------------\n');
fprintf(1,'Beginning %s()\n',subroutine_name);
fprintf(1,'-------------------------------------------------------------\n');

tmp  = [res.par];
TE   = [tmp.TE];
nt   = numel(TE);
if (numel(unique(TE)) < (nparams+1)), res = []; fprintf(2,'ERROR: T2 fit needs at least %1d different TE values\n',(nparams+1)); return; end
amps     = [res(:).amp];
amps     = reshape(amps, peaks.count,nt);   % pull our peak fit amps from result struct array
ampfit   = amps * 0;
ampfixed = amps * 0;
TEfine   = linspace(0,max(TE(:)),nt*10);
finefit  = zeros(peaks.count,nt*10);
for i = 1:peaks.count
    [~,fit,fitparams,rnorm,fitfine,datafixed] = t2fit(TE,amps(i,:),TEfine,nparams);
    res(1).t2fit_M0(i)      = fitparams(1);
    res(1).t2fit_t2(i)      = fitparams(2);
    res(1).t2fit_offset(i)  = fitparams(3);
    res(1).t2fit_rnorm(i)   = rnorm;
    res(1).t2fit_rsq(i)     = corr2(amps(i,:),fit)^2;
    ampfit(i,:)             = fit;
    ampfixed(i,:)           = datafixed;
    finefit(i,:)            = fitfine;
end
% return the fitted amps and TEs
for j = 1:nt
    res(j).TE       = TE(j);   
    res(j).ampfit   = ampfit(:,j)';
    res(j).ampfixed = ampfixed(:,j)';    % this is for compatibility w/ irfit() results
end
res(1).TEfine  = TEfine;
res(1).finefit = finefit;  % this is a bit of a kluge
res(1).fittype = 't2decay';
end


