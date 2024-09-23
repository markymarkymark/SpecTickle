function [fit,Fits,Pars,z,a] = me_hsvd(fid,nt,nz,bw,skip,verbose,do_hlsvd)
%HSVD State-space method for the estimation of an FID signal.
%
%  [FIT,FITS,PARS,Z,A] = HSVD(FID,NT,NZ,BW) returns the best estimate
%  to an FID signal using the the state-space singular value
%  decomposition algorithm.
%
%  FID is the data vector. NT is the number of points in fid to use
%  for fitting. this can be <= length(fid), but fit will always be
%  the same length as fid. NZ is the number of signal poles to find,
%  although the actual number of signal poles returned (see z in
%  "outputs" below) may be less if there are poles found outside the
%  unit circle (ie, exponentially increasing signals are discarded).
%  BW is the bandwidth of fid in Hz.
%  
%  FIT is the best estimate to FID. FITS is an NxM matrix of
%  component fits, where N is the length of fit and M is the length
%  of z (see below). Z is vector of signal poles and A is vector of
%  signal amplitudes. PARS is an Mx4 matrix of model function
%  parameters, where M is the length of z (see below). The 4
%  parameters are stored as follows:
%
%        PARS(:,1) = frequencies in Hz
%        PARS(:,2) = amplitudes in arbitrary units
%        PARS(:,3) = T2s in msec
%        PARS(:,4) = phase in degrees
%
%  REFERENCE: H Barkhuijsen et al, J Magn Reson 73:553-557 (1987).
%

if (nargin < 7), do_hlsvd = 0; end   % HLSVD
if (nargin < 6), verbose  = 0; end   % print stuff
if (nargin < 5), skip     = 0; end   % ignore some initial FID points
if (skip > 0)
    if (verbose), fprintf(1,'    Ignoring %1d initial FID points for HSVD decomp\n',skip); end
    np   = size(fid,1);
    Xfid = fid(skip+1:np);
else
    Xfid = fid;
end

% initialize a nearly square Hankel matrix
m = round(nt/2);
l = nt-m+1;
H = zeros(l,m);
if (verbose), fprintf(1,'    Hankel matrix is %1dx%1d\n',l,m); end

% fill the normalized Hankel data matrix (column-by-column)
x = Xfid/max(abs(Xfid));
for k=1:m, H(:,k) = x(k:k+l-1); end

% SVD the Hankel data matrix
if (~do_hlsvd)
    [U,S,V] = svd(H,0);
else
    [U,S,V] = svds(H,l);  % maybe this is Lanczos?
end

% form the top and bottom truncated left singular vector matrices
Ut = U(2:l,1:nz);
Ub = U(1:l-1,1:nz);

% solve for Z' by least squares
Zp = pinv(Ub) * Ut;

% find the eigenvalues of Z'
z = eig(Zp);

% remove poles outside the unit circle (ie, exponentially increasing signals)
z = z(abs(z) < 1);
nz = size(z,1);

% sort poles by increasing frequency
[ord,k] = sort(angle(z));
z = z(k);

% form the Z matrix and s data vector
s = Xfid(1:nt);
Z = zeros(nt,nz);
t = (0:nt-1)' + skip;
for k=1:nz, Z(:,k) = z(k).^t; end

% find the amplitudes by least squares
a = pinv(Z) * s;
%a = Z \ s;

% generate the composite fit and component fits
nf = size(fid,1);
t = (0:nf-1)';
Fits = zeros(nf,nz);
for k=1:nz, Fits(:,k) = a(k)*z(k).^t; end
fit = sum(Fits,2);

% calculate the canonical parameters: freq, t2, amp, and phase
freq = bw*angle(z)/(2*pi);
t2 = -1e3./(log(abs(z))*bw);
amp = abs(Fits(1,:))';
phase = angle(Fits(1,:))'.*(180/pi);
Pars = [freq, amp, t2, phase];

% show the results
if (verbose && 0)
    disp('         freq(Hz)    amp(au)    T2*(ms)   phi(deg)');
    for k=1:nz, fprintf(1,'%2d) %10.2f %10.2f %10.2f %10.2f\n', k,Pars(k,:)); end
end
