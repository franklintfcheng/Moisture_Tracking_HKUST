function poly = draw_polyshape(varargin)
        
    pnames = {'type', 'center', 'radius'};

    dflts  = cell(1, length(pnames));

              [type,   center,   radius]...
                           = internal.stats.parseArgs(pnames, dflts, varargin{:});

    %%
    n = 100;
    if isequal(type, 'circle')
        theta = (0:n-1)*(2*pi/n);
        x = center(1) + radius*cos(theta);
        y = center(2) + radius*sin(theta);
        poly = polyshape(x,y);
    else
        error('code not set for %s type of polyshape!', type)
    end
end