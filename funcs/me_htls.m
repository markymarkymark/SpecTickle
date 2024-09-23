function [fit,Fits,Pars,q,a] = me_htls(fid,nt,nz,bw,skip,verbose)

%function [frequencies, dampings, basis, ahat] = htls(y, fs)

%
% Decompose the signal fid using the method of Van Huffel, et al. (1993)
%
% Obligatory arguments:
% 'fid' is the FID (a linear combination of complex exponentials).
% 'bw' is the sampling frequency (bandwidth).
%
% Outputs:
% 'frequencies' - frequencies of components in signal
% 'dampings' - damping of components in signal
% 'basis' - basis vectors, one for each component
% 'ahat' - amplitudes of each basis in signal
%
% Author: Greg Reynolds (remove.this.gmr001@bham.ac.uk)
% Date: August 2006 

if (nargin < 6), verbose = 0; end   % print stuff
if (nargin < 5), skip = 0; end      % ignore some initial FID points
if (skip > 0)
    if (verbose), fprintf(1,'    Ignoring %1d initial FID points for HSVD decomp\n',skip); end
    np   = size(fid,1);
    Xfid = fid(skip+1:np);
else
    Xfid = fid;
end


%%N = length(fid); % original code
N = nt;
L = floor(0.5*N);
M = N+1-L;

% H is the LxM Hankel LP matrix
H = hankel(Xfid(1:L), Xfid(L:N));

[U,S,V] = svd(H);
s = diag(S);

% it is always better to overestimate apparently
%% = estimate_model_order(s, N, L)+10;   % original code
K = nz;

% construct Uk
Uk = U(:,1:K);

% find Ukt and Ukb
Ukt = discardrow(1, Uk);
Ukb = discardrow(size(Uk,1),Uk);

[U,S,V] = svd([Ukb Ukt]);
V12 = V(1:K, K+1:2*K);
V22 = V(K+1:2*K, K+1:2*K);

Zp = -V12*inv(V22);

% find the modes
q = eig(Zp);
q = log(q);

dt = 1 /bw;
freq = imag(q)/(2*pi) / dt;
dampings = real(q) / dt;

% sort poles by increasing frequency
[freq,ifreq] = sort(freq);
dampings = dampings(ifreq);

% construct the basis
nf = size(fid,1);
t = (0:dt:(nf-1)*dt);
basis = exp(t.'*(dampings.'+1i*2*pi*freq.'));

% compute the amplitude estimates
a = pinv(basis)*fid;

% generate the composite fit and component fits
t = (0:nf-1)';
Fits = zeros(nf,nz);
for k=1:nz, Fits(:,k) = a(k)*basis(:,k); end
fit = sum(Fits,2);

% calculate the canonical parameters: freq, t2, amp, and phase
t2 = -1e3./dampings;
amp = abs(Fits(1,:))';
phase = angle(Fits(1,:))'.*(180/pi);
Pars = [freq, amp, t2, phase];

end
