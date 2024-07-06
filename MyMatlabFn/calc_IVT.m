function [uIVT, vIVT] = calc_IVT(varargin)

    %/ IMPORTANT: data has to be 4D, with the 3rd dim to be the level dim.
    pnames = {'U', 'V', 'q', 'level', 'level_unit', 'lnP_coor', 'topo_lon', 'topo_lat', 'sp_mode', 'sp', 'interp_mode'}; 
    dflts  = repmat([], length(pnames), 1);
    [          U,   V,   q,   level,   level_unit,   lnP_coor,   topo_lon,   topo_lat,   sp_mode,   sp,   interp_mode] = internal.stats.parseArgs(pnames, dflts, varargin{:});
%%
    %==========================================================================
    %/      Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last Update: 17 May 2024
    %/
    %/ Description: This function is designd to calculate vertically
    %/              integrated vapor transport (IVT)
    %==========================================================================
    
    qu = q.*U;
    qv = q.*V;
    wd = 997;  %/ density of water 
    
    uIVT = 1000/wd*vertinte('data', qu, 'level', level, 'level_unit', level_unit, 'lnP_coor', lnP_coor, 'topo_lon', topo_lon, 'topo_lat', topo_lat, 'sp_mode', sp_mode, 'sp', sp, 'interp_mode', interp_mode);
    vIVT = 1000/wd*vertinte('data', qv, 'level', level, 'level_unit', level_unit, 'lnP_coor', lnP_coor, 'topo_lon', topo_lon, 'topo_lat', topo_lat, 'sp_mode', sp_mode, 'sp', sp, 'interp_mode', interp_mode);

end