
function sequence = sequence(varargin)

    pnames = {'data','skipgap'};
    dflts  = {[], []};
    [data, skipgap] = internal.stats.parseArgs(pnames, dflts, varargin{:});
    
    % [n1, n2] = size(data);
    
    
    if size(data, 1) == 1   data = data';  end  %/ convert to column vector if it's not so
    
    if any(ismember(diff(data), -1))
        error('The input vector is not strictly increasing. Please check.');
    end
    
    if isempty(skipgap)
        nodes = ~ismember(diff(data), 1); %/ discontinued nodes
        
    else
        nodes = ~ismember(diff(data), 1:skipgap+1); %/ discontinued nodes + skip small gaps
    %     nodes = diff(data)~= [1, skipgap+1];  
        disp('NOTE: the skipped days will not be filled into the input vector, just to merge sequences.');
    end
    firstInd = [true; nodes];
    
    sequence = {};
    e = 0;
    for i = 1:length(data)
        if firstInd(i) == 1  %/ if reach the discontinued point, create a new cell
            e = e + 1; 
            sequence{e, 1} = [];
        end
        sequence{e, 1} = [sequence{e}, data(i)];
    end

end