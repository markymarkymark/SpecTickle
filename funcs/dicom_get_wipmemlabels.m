function [WipStruct,LongLabels,DoubleLabels] = dicom_get_wipmemlabels(longvals,doublevals,force_double)

if (nargin < 3), force_double = 1; end
WipStruct = []; LongLabels = []; DoubleLabels = [];

% if cell arrays, then NaNs mean no value was in header - assume = 0!
if (iscell(longvals))
    longvals(cellfun(@isnan,longvals))     = {int32(0)};
    doublevals(cellfun(@isnan,doublevals)) = {0};
    longvals   = cell2mat(longvals);
    doublevals = cell2mat(doublevals);
end
max_longs     = 20;
max_doubles   = 16;
LabelID_index = max_longs + max_longs + max_doubles + 1;

% --- Did we encode param names in the unused LONG values? ---
if (numel(longvals) >= LabelID_index && ~isnan(longvals(LabelID_index)))
    label = char(typecast(longvals(LabelID_index),'uint8'));
    if (isequal(label,'PENN'))
        LongLabels   = cell(max_longs,1);
        DoubleLabels = cell(max_doubles,1);
        for i = 1:max_longs
            LongLabels{i} = char(typecast(longvals(i+max_longs),'uint8'));
        end
        for i = 1:max_doubles
            DoubleLabels{i} = char(typecast(longvals(i+2*max_longs),'uint8'));
        end
        for i = 1:max_longs
            if (~isequal(LongLabels{i},'xxxx'))
                if (~force_double), WipStruct.(LongLabels{i}) = longvals(i); % these are int32, which MATLAB complains about alot
                else,               WipStruct.(LongLabels{i}) = double(longvals(i)); end
            end
        end
        for i = 1:max_doubles
            if (~isequal(DoubleLabels{i},'xxxx'))
                WipStruct.(DoubleLabels{i}) = doublevals(i);
            end
        end
    end
end
end