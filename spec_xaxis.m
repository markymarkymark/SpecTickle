function [t,hz,ppm] = spec_xaxis(np,bw,f0,tshift,ppmshift)
% spec_xaxis - calculate time and frequency x-axis vectors for spectroscopy
% 
% Syntax:
%   [t,hz,ppm] = spec_xaxis(np,bw,f0,[tshift,ppmshift])
%
% Created: Mark A. Elliott, PhD
%   melliott@upenn.edu
%   https://www.med.upenn.edu/CAMIPM/mark-elliott.html
%
if (nargin < 4 || isempty(tshift)),   tshift   = 0; end
if (nargin < 5 || isempty(ppmshift)), ppmshift = 0; end

t   = ((0:np-1)')/bw + tshift;
hz  = (0:np-1)';
hz  = (hz/np-0.5)*bw;
ppm = -hz/f0 + ppmshift;
end
