function [t2,fit,pars,rnorm,fitfine,fitdata] = t2fit(te,data,tefine,nparams)
% Fit multi-TE data for T2
% Can fit either 2 or 3 param model. 3rd param is saturation efficiency
%
% Syntax:
%   [t2,fit,pars,rnorm,fitfine] = t2fit(ts,data,[tefine],[nparams])
%
% Created: Mark A. Elliott, PhD
%   melliott@upenn.edu
%   https://www.med.upenn.edu/CAMIPM/mark-elliott.html

if (nargin < 3), tefine = []; end
if (nargin < 4), nparams = 2; end
opts = optimoptions(@lsqcurvefit,'Display','off');

if (nparams == 3)
    pguess = [data(1), 20, data(end)];
    [pars,rnorm] = lsqcurvefit(@t2func_3param,pguess,te,data,[],[],opts);
    fit = t2func_3param(pars,te);
    if (~isempty(tefine)), fitfine = t2func_3param(pars,tefine);  end
else
    pguess = [data(1), 20];
    [pars,rnorm] = lsqcurvefit(@t2func_2param,pguess,te,data,[],[],opts);
    fit = t2func_2param(pars,te);
    if (~isempty(tefine)), fitfine = t2func_2param(pars,tefine);  end
    pars(3) = 0.0;  % force assumed DC offset of 0
end
t2 = pars(2);
fitdata = data; % this is for compatibility w/ irfit(), in case the fitted data was changed (e.g. sorted)
end

% ------------------------------------------------------------------------
function F = t2func_2param(p,x)
F = p(1)*exp(-x/p(2));
end

% ------------------------------------------------------------------------
function F = t2func_3param(p,x)
F = p(1)*exp(-x/p(2)) + p(3);
end