function spec = m1dfft(fid,reverse)
% perform 1D FFT along FIRST dimension of FID

if (nargin < 2), reverse = 0; end
[np,nfids] = size(fid);
spec       = complex(zeros(np,nfids),zeros(np,nfids));

for i=1:nfids
    if (~reverse)
        fid(1,i)  = fid(1,i)/2;
        spec(:,i) = fftshift(fft(fid(:,i)));
    else
        spec(:,i) = ifft(ifftshift(fid(:,i)));
        spec(1,i) = spec(1,i)*2;
    end
end

end

