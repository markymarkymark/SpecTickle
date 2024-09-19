function val = checkstruct(s,fieldname,defval,fieldname2)

if (nargin == 2) % this is called like checkstruct(struct1,struct2)
    s2 = fieldname;
    names = fieldnames(s2);
    for i=1:numel(names)
        if (~isfield(s,names{i})), s = setfield(s,names{i},getfield(s2,names{i})); end
    end
    val = s;
    return
end

if (~isfield(s,fieldname))
    val = defval;
    return
end

if (nargin > 3)
    s2 = getfield(s,fieldname);
    if (~isfield(s2,fieldname2))
        val = defval;
        return
    else
        s = s2;
        fieldname = fieldname2;
    end
end

val = getfield(s,fieldname);
%%if (ischar(val) && isempty(val)), val = defval; end  % empty strings get the default value
if (isempty(val)), val = defval; end  % empty vars get the default value

return
