%%
function array_split = job_distributor(varargin)

    pnames = {'array', 'Njob', 'job_id'};  
    dflts  = cell(1, length(pnames));

    [          array,    Njob,   job_id] ...
                   = internal.stats.parseArgs(pnames, dflts, varargin{:}); %/ parse function arguments

    %%
    if ~isvector(array)  error('The input array must be an vector only!');  end
    if job_id > Njob     error('job_id must not be larger than Njob!');     end
    
    %/ split the job
    f = floor(length(array)/Njob);

    if job_id == Njob
        array_split = (1 + f*(job_id-1)):length(array); %/ to workaround the uneven job split.
    else
        array_split = (1:f) + f*(job_id-1);
    end
            

end