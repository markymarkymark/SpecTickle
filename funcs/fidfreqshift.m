function [fshift,fit,pars,rnorm] = fidfreqshift(t,fid,fidref,pguess,Nfit,skippts)
% Find freq shift in Hz to match FID to FIDREF
% Can fit one param for freq shift, or 2 for shift and amp

global fidx   % this is the FID to be shifted, need to share with optim. func

if (nargin < 5), Nfit    = size(fid,1); end % fit to only the first Nfit points if desired
if (nargin < 6), skippts = 0; end           % ignore initial FID points

if (~isequal(class(fid),   'double')), fid    = double(fid);    end
if (~isequal(class(fidref),'double')), fidref = double(fidref); end

fidx       = fid((1+skippts):(Nfit+skippts));                                                       % need to share w/ optimization func
fidref_cat = [real(fidref((1+skippts):(Nfit+skippts))) ; imag(fidref((1+skippts):(Nfit+skippts)))]; % cat means convert complex to concatened real/imag
time       = t((1+skippts):(Nfit+skippts));

% --- Do the fit ---
opts = optimoptions(@lsqcurvefit,'Display','off');
switch numel(pguess)
    case 1   % fit freq shift
        func = @fidfreqshift_func1;
    case 2   % fit freq shift and amp
        func = @fidfreqshift_func2;
    case 3   % fit freq shift, amp & phase
        func = @fidfreqshift_func3;
end
[pars,rnorm] = lsqcurvefit(func,pguess,time,fidref_cat,[],[],opts);

% --- Get the results ---
fshift  = pars(1);
fidx    = fid;      % if fit was done on truncated FID, return fit to full FID
fid_cat = func(pars,t); 
np      = numel(fid_cat);
fit     = complex(fid_cat(1:np/2),fid_cat((np/2+1):end)); % turn concatenated real result into complex
end

% ------------------------------------------------------------------------
function fid_cat = fidfreqshift_func1(p,t)
global fidx
fidguess  = fidx .* exp(-1i*p(1)*2*pi*t);
fid_cat   = [real(fidguess)  ; imag(fidguess)];
end

% ------------------------------------------------------------------------
function fid_cat = fidfreqshift_func2(p,t)
global fidx
fidguess  = p(2) * fidx .* exp(-1i*p(1)*2*pi*t);
fid_cat   = [real(fidguess)  ; imag(fidguess)];
end

% ------------------------------------------------------------------------
function fid_cat = fidfreqshift_func3(p,t)
global fidx
fidguess  = p(2) .* exp(-1i*p(3)) .* exp(-1i*p(1)*2*pi*t) .* fidx;
fid_cat   = [real(fidguess)  ; imag(fidguess)];
end

