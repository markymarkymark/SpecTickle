function [fwhm,s,s_interp,hz_interp,pfwhm1,pfwhm2] = fwhm(spec,hz,interp_factor,do_0phase)
% FWHM - Caculate the FWHM of a spectrum with a single peak in it.
% Uses interpolation to get a better resolution estimate than the sampling
% 
% Syntax:
%   [fwhm,s,s_interp,hz_interp,pfwhm1,pfwhm2] = fwhm(spec,[hz,interp_factor,do_0phase])
% Inputs:
%   spec          - real or complex spectrum with single peak in it
%   hz            - x-axis in Hz
%   interp_factor - factor to interpolate for higher resolution (default 20)
%   do_0phase     - first auto 0-order phase the (complex!) spectrum (default = 1)
%
% Outputs:
%   fwhm       - the FWHM value (in Hz)
%   s         - the (possibly phased) real spectrum
%   s_interp  - the higer spectral resolution interpolated spectrum
%   hz_interp - the highe resolution x-axis
%
% Created: Mark A. Elliott, PhD
%   melliott@upenn.edu
%   https://www.med.upenn.edu/CAMIPM/mark-elliott.html
%

if (nargin < 2 || isempty(hz)),            hz            = 1:size(spec,1); end
if (nargin < 3 || isempty(interp_factor)), interp_factor = 20;  end
if (nargin < 4 || isempty(do_0phase)),     do_0phase     = 1;  end

hz_step    = hz(2)-hz(1);         % need to increase spectral resolution to get better FWHM estimate
hz_interp  = hz(1):hz_step/interp_factor:hz(end);
if (do_0phase), s = phase_zero(spec); end
s          = real(s);
s          = s/max(s(:));
s_interp   = spline(hz,s,hz_interp);
[~,pmax]   = max(s_interp);
[~,pfwhm1] = min(abs(s_interp(1:pmax) - 0.5));
[~,pfwhm2] = min(abs(s_interp(pmax:end) - 0.5));
pfwhm2     = pmax + pfwhm2 - 1;
fwhm       = hz_interp(pfwhm2)-hz_interp(pfwhm1);
end