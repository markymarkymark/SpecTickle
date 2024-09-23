function [t1,fit,pars,rnorm,fitfine,fitdata] = irfit(ti,data,tifine,nparams)
% Fit Inversion Recovery data for T1
% Can fit either 2 or 3 param model. 3rd param is inversion efficiency
% This function negates data points prior to min(data), and tries 2 fits, with +/_ data(imin)
%
% Syntax:
%   [t1,fit,pars,rnorm,fitfine,fitdata] = irfit(ts,data,[tifine],[nparams])
%
% Created: Mark A. Elliott, PhD
%   melliott@upenn.edu
%   https://www.med.upenn.edu/CAMIPM/mark-elliott.html

if (nargin < 3), fitfine = []; end
if (nargin < 4), nparams = 2;  end
opts = optimoptions(@lsqcurvefit,'Display','off');

% sort for increasing TI
[ti,iTI] = sort(ti);
data     = data(iTI);

% make 2 data sets, with min(abs(data)) being pos or neg
[~,imin] = min(data(:));
if (imin > 1)
    data(1:imin-1) = -data(1:imin-1);  % make points before zero-crossing negative
end
data2       = data;
data2(imin) = -data2(imin);          % try 2nd version of data w/ sign of zero-crossed negated

% --- 2 parameter IR fit ---
% pguess = [1000, data(end)];
% [pars2,rnorm] = lsqcurvefit(@irfunc_2param,pguess,ti,data,[],[],opts);
% t1 = pars2(1);
% fit2 = t1func_2param(pars2,ti);

if (nparams == 3)
    pguess = [1000, data(end), 1];
    [pars,rnorm] = lsqcurvefit(@irfunc_3param,pguess,ti,data,[],[],opts);
    fit = irfunc_3param(pars,ti);

    % try fit where zero cross point is sign-reversed
    [parsx,rnormx] = lsqcurvefit(@irfunc_3param,pguess,ti,data2,[],[],opts);
    fitx = irfunc_3param(parsx,ti);
end

% use fit to data2 if better residue
if (rnormx < rnorm)
    fit     = fitx;
    pars    = parsx;
    rnorm   = rnormx;
    fitdata = data2;
else
    fitdata = data;
end
t1 = pars(1);
if (nargin > 2), fitfine = irfunc_3param(pars,tifine); end
% undo the sort
fit(iTI)     = fit;
fitdata(iTI) = fitdata;
end


% ------------------------------------------------------------------------
% f = M0*(1-2*exp)
function F = irfunc_2param(p,x)
F = p(2)*(1 - 2*exp(-x/p(1)));
end

% ------------------------------------------------------------------------
% f = M0*(1-2*k*exp)
function F = irfunc_3param(p,x)
F = p(2)*(1 - 2*p(3)*exp(-x/p(1)));
end