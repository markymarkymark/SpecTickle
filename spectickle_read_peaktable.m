function [d,peakfile] = spectickle_read_peaktable(peakfile)
% Reads the tab-delimited table of resonances, fieldstrength, nucleus, etc, ...

d = [];
peak_table_default = 'spectickle_peaktable.tsv';

% --- if no file provided, get default peak file ---
if (nargin < 1),        peakfile = ''; end
if (isempty(peakfile))
    thisfile = mfilename('fullpath');
    peakfile = fullfile(fileparts(thisfile), 'tables', peak_table_default);
end
if (~isfile(peakfile)), return; end % error
fprintf(1,'Loading peak definitions from %s\n',peakfile);

% --- Get structure of table file ---
d = tdfread(peakfile);

% --- Remove blank rows from all fields ---
notblank = ~isnan(d.ppm);   % find blank rows
f        = fieldnames(d);
nf       = numel(f);
for i=1:nf
    vals     = d.(f{i});
    if (ischar(vals)), vals = deblank(string(vals)); end
    d.(f{i}) = vals(notblank);
end
nv = numel(d.FieldStrength);
end






