%%

function [c, s, double_check] = trig_encoding(varargin)
    
    pnames = {'x', 'T'};
    dflts  = cell(1, length(pnames));
    [          x,   T] =  internal.stats.parseArgs(pnames, dflts, varargin{:});
    
    %/ This function is to transform the periodic attribute (e.g., 24 hours of
    %/   day, directional heading, etc.) into two encoded variables to
    %/   adequately describe the distance in the values.
    %/
    %/ It is used prior to distance-based clustering like the K-means.
    
    %/ Its disadvantages were discussed in Vejmelka et al (2009).
    
    %/ x = the periodic variable
    %/ T = the period of x
    
    if size(x, 1) ~= 1 && size(x, 2) ~= 1   error('This function only handles vectors!');   end
    if size(x, 2) ~= 1    x = x';   end
    
    c = cos(2*pi*x./T); %/ mapping
    s = sin(2*pi*x./T); %/ mapping
    
    %/ Output a matrix for a double-check.
    double_check = [x, c, s];
    double_check = round(double_check, 4); 


end