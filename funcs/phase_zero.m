function [pspec,phi,pspecx] = phase_zero(spec,p1,p2,xvec)
%  Auto zero-order phase a 1D spectrum
%
% Syntax:
%  [pspec,phi,pspecx] = phase_zero(spec,[p1,p2],[xvec])  
%
% Created: Mark A. Elliott, PhD
%   melliott@upenn.edu
%   https://www.med.upenn.edu/CAMIPM/mark-elliott.html
pspec = []; phi = []; pspecx = [];
switch (nargin)
    case 1
        specx = spec;
    case 3
        specx = spec(p1:p2);
    case 4        
        [~,i1] = min(abs(xvec(:) - p1));
        [~,i2] = min(abs(xvec(:) - p2));
        specx = spec(i1:i2);
    otherwise
        fprintf(2,'ERROR: %s() - wrong number of arguments\n',mfilename());
        return
end

phix  = 0:1:360;
pi    = 3.14159;
smax  = max(abs(specx(:)));
specx = specx/smax;
nphix = numel(phix);
maxarea = -1e9;
for i=1:nphix
    pspec = specx * exp(+1i*phix(i)*pi/180);
%    area  = sum(real(pspec(:)));
    area  = max(real(pspec(:)));
    if (area > maxarea)
        phi = phix(i);
        maxarea = area;
    end
end
pspecx = smax * specx * exp(complex(0,phi*pi/180)); 
pspec  =         spec * exp(complex(0,phi*pi/180));     % apply phase cor to whole input spectrum
end

