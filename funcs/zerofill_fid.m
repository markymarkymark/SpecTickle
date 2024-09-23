function [fidfilled,npout,npin,nfids] = zerofill_fid(fid,npfilled)

fidfilled = [];
if (nargin < 2), return; end

[npin,nfids] = size(fid);
if (npfilled > 0)  % positive number means fill to this size
    if (npfilled <= npin), fprintf(2,'ERROR: requested zero filled points (%1d) is less than FID size (%1d)\n',npfilled,npin); return; end
    zfill = zeros(npfilled-npin,nfids);  
else                % negaive number means fill by this |factor|
    zfill = zeros((abs(npfilled)-1)*npin,nfids);  
end
fidfilled = [fid ; zfill];
npout     = size(fidfilled,1);
end
