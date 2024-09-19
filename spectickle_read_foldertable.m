function [d,folderfile] = spectickle_read_foldertable(folderfile)
% Reads the tab-delimited table folder matching strings and corresponding experiemnt type

d = [];
folder_table_default = 'spectickle_foldertable.tsv';

% --- if no file provided, get default peak file ---
if (nargin < 1), folderfile = ''; end
if (isempty(folderfile)), folderfile = which(folder_table_default); end
if (isempty(folderfile)), folderfile = folder_table_default; return; end % error
fprintf(1,'Loading folder matching definitions from %s\n',folderfile);

% --- Get structure of table file ---
d = tdfread(folderfile);

% --- convert strings to structure of cell strings ---
f        = fieldnames(d);
nf       = numel(f);
for i=1:nf
    vals     = d.(f{i});
    vals     = deblank(string(vals));
    d.(f{i}) = cellstr(vals);
end
nv = numel(d.Folder_Match);
end






