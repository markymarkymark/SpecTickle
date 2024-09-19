function [t1,fit,pars,rnorm,fitfine,fitdata] = srfit(ts,data,tsfine,nparams)
% Fit Saturation Recovery or Progressive Saturation data for T1
% Can fit either 2 or 3 param model. 3rd param is saturation efficiency
%
% Syntax:
%   [t1,fit,pars,rnorm,fitfine] = srfit(ts,data,[tsfine],[nparams])
%
% Created: Mark A. Elliott, PhD
%   melliott@upenn.edu
%   https://www.med.upenn.edu/CAMIPM/mark-elliott.html

if (nargin < 3), tsfine = []; end
if (nargin < 4), nparams = 3; end
opts = optimoptions(@lsqcurvefit,'Display','off');

% --- 3 param model ---
if (nparams == 3)
    pguess = [1000, data(end), 1.0];                   
    [pars,rnorm] = lsqcurvefit(@srfunc_3param,pguess,ts,data,[],[],opts);
    fit = srfunc_3param(pars,ts);
    if (~isempty(tsfine)), fitfine = srfunc_3param(pars,tsfine); end
% --- 2 param model ---
else
    pguess = [1000, data(end)];
    [pars,rnorm] = lsqcurvefit(@srfunc_2param,pguess,ts,data,[],[],opts);
    fit = srfunc_2param(pars,ts);
    if (~isempty(tsfine)), fitfine = srfunc_2param(pars,tsfine); end
    pars(3) = 1.0;  % force assumed SR efficiency of 1.
end
t1 = pars(1);
fitdata = data; % this is a place holder, see irfit()
end

% ------------------------------------------------------------------------
function F = srfunc_2param(p,x)
F = p(2)*(1 - exp(-x/p(1)));
end

% ------------------------------------------------------------------------
function F = srfunc_3param(p,x)
%F = p(2) - p(3)*exp(-x/p(1));
F = p(2)*(1 - p(3)*exp(-x/p(1)));
end