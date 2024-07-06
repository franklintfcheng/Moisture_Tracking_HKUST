function seg = get_segment(varargin)
    
    %/ create a set of valid parameters and their default value
    pnames = {'logical_array', 'n'};  
    dflts  = cell(1, length(pnames));
    
    [          logical_array,   n] ...
                   = internal.stats.parseArgs(pnames, dflts, varargin{:}); %/ parse function arguments

%%
    %/===========================================================
    %/ Authoer: Franklin Cheng
    %/ Last update: Jan 19, 2024
    %/
    %/ NOTE: This function is designed to get segments with length >= n 
    %/       from a logical array.
    %/    
    %/===========================================================

    %/ TESTING
    % logical_array = [0 1 1 1 0 1 0 1 0 1 0]';
    % logical_array = [1 1 1 1 0 0 0 1 1 1 0]';
    % logical_array = [1 1 1 1 1 0 1 1 1 1 1]';
    % logical_array = [0 1 1 1 0 0 0 1 1 1 1]';

    %/ Output all segments with the minimum length of 1
    if isempty(n)  
        n = 1;   
    end

    %/ Make it a column vector
    if size(logical_array, 1) == 1  
        logical_array = logical_array';  
    end  

    %/ Append a fake, opposite index at the beginning before diff()
    q = logical_array;
    if logical_array(1) == 1
        q = [0; q];
    else
        q = [1; q];
    end

    d      = diff(q);
    st_ind = find(d ==  1); 
    ed_ind = find(d == -1) - 1;  %/ Shift 1 index forward to get the correct end point
    
    if logical_array(1) == 0
        ed_ind(1) = [];
    end
    if logical_array(end) == 1
        ed_ind = [ed_ind; length(logical_array)];
    end

    seg_gte_n   = find((ed_ind-st_ind+1) >= n);
    seg         = [st_ind(seg_gte_n), ed_ind(seg_gte_n)];  %/ Output seg with (starting index, ending index)
end