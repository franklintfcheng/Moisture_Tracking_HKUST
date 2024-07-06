%%
function vertices_new = close_vertices(varargin)
    %/ Author: Fandy
    %/ Last Update: 28 May, 2022
    
    switch nargin
        case 1  %/ given 1 input
            vertices = varargin{1};
        otherwise
            error('Unexpected inputs')
    end

    %/ For coding convenience, first append vertices with nans if not so
    if ~all(isnan(vertices(end,:)))
        vertices = [vertices; nan nan];
    end
    
    %/ Find the separation point (by nan)
    ind_nan = find(isnan(vertices(:,1)));
    
    %/ Get the start and end point of each segment
    ind_st = [1; ind_nan(1:end-1)+1];
    ind_ed = ind_nan;
    
    %/ Insert the start point at the end of each segment -> closing the vertices
    %/ Append to vertices_new 
    vertices_new = [];
    for i = 1:length(ind_st)
        seg = vertices(ind_st(i):ind_ed(i)-1,:);
        seg = [seg; seg(1,:); [nan, nan]];
        vertices_new = [vertices_new; seg];
    end
    
end