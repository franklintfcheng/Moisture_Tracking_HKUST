%%
function rename_files(varargin)

%/ create a set of valid parameters and their default value
pnames = {'filepath', 'keyword', 'keyword_new'};  
dflts  = cell(1, length(pnames));

[          filepath, keyword, keyword_new] ...
               = internal.stats.parseArgs(pnames, dflts, varargin{:}); %/ parse function arguments

fstruct    = dir([filepath, '*.*']);  %/ since it auto sorts in ascending order, so no need to do anything.
file_list  = {fstruct.name}';

a = strfind(file_list, keyword);
ind = find(~cellfun(@isempty,a));

if isempty(ind)   warning('None of the files contains ''%s'' to be replaced with!', keyword);   end

for i = 1:length(ind)
    fprintf('*** %d/%d ***\n', i, length(ind))
    oldname = file_list{ind(i)};

    st_ind_char = a{ind(i)};

    prefix = oldname(1:st_ind_char-1);                          %/ char string before the keyword.
    suffix = oldname(st_ind_char+length(keyword):end);        %/ char string after the keyword.

    newname = [prefix, keyword_new, suffix];

    if isfile([filepath, oldname])
        movefile([filepath, oldname], [filepath, newname]);
        fprintf('!!! Renamed from %s to %s !!!\n', oldname, newname);
    else
        error('expected file "%s" not found\n', oldname);
    end
end


end