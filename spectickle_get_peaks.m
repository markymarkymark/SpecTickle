function peaks = spectickle_get_peaks(peaktable,opt)

peaks = [];

% --- Find all defined peak groups for given field, tissue and nucleus ---
%p = (strcmpi(peaktable.FieldStrength,opt.fieldstrength) & strcmpi(peaktable.Nucleus,opt.nucleus) & strcmpi(peaktable.Tissue,opt.tissue_type));
p = (contains(peaktable.FieldStrength,opt.fieldstrength,'IgnoreCase',true) & contains(peaktable.Nucleus,opt.nucleus,'IgnoreCase',true) & contains(peaktable.Tissue,opt.tissue_type,'IgnoreCase',true));
if (~any(p)), fprintf(2,'ERROR: Did not find any defined peak groups for FieldStrength/Tissue/Nucleus = %s/%s/%s\n',opt.fieldstrength,opt.tissue_type,opt.nucleus); return; end

% --- Now find MATCHED peak groups for given field, tissue and nucleus ---
nv = numel(p);
p = zeros(nv,1,'logical');
for i=1:numel(opt.peak_groups)
    px = (contains(peaktable.FieldStrength,opt.fieldstrength,'IgnoreCase',true) & contains(peaktable.Nucleus,opt.nucleus,'IgnoreCase',true) & contains(peaktable.Tissue,opt.tissue_type,'IgnoreCase',true) & strcmpi(peaktable.Group,opt.peak_groups{i}));
    if (~any(px)), fprintf(2,'ERROR: Did not find peak group "%s" for FieldStrength/Tissue/Nucleus = %s/%s/%s\n',opt.peak_groups{i},opt.fieldstrength,opt.tissue_type,opt.nucleus); return; end
    p = p | px;
end

% --- Fill peak structure for fitting ---
peaks.labels    = cellstr(peaktable.Label(p));
peaks.groups    = cellstr(peaktable.Group(p));
peaks.ppm       = peaktable.ppm(p);
peaks.ppm_tol   = peaktable.ppmtol(p);
peaks.lw_min    = peaktable.linewidth_min(p);
peaks.lw_max    = peaktable.linewidth_max(p);
peaks.phase_min = -peaktable.phase_tol(p);
peaks.phase_max = peaktable.phase_tol(p);
peaks.yoke_group= cellstr(peaktable.yoke_group(p));
peaks.yoke_amp  = peaktable.yoke_amp(p);
peaks.count     = sum(p);

if (opt.peaktable_shift ~= 0)
    if (opt.verbose), fprintf(1,'Shifting frequency locations of all peaks to be fit by %1.2f ppm\n',opt.peaktable_shift); end
    peaks.ppm = peaks.ppm + opt.peaktable_shift;
end
end
