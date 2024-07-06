%%

function [bndry_data] = point2bndry(varargin)

    pnames = {'point_lon', 'point_lat', 'res_lon', 'res_lat'};
    dflts  = cell(length(pnames), 1);
    [           point_lon,   point_lat,   res_lon,   res_lat] = internal.stats.parseArgs(pnames, dflts, varargin{:});

    if isempty(res_lon)  error('Please define the resolution of lon!'); end
    if isempty(res_lat)  error('Please define the resolution of lat!'); end
    
    bndry_data = cell(length(point_lon), 1);
    for i = 1:length(point_lon)
        bndry_data{i} = [point_lon(i) - res_lon/2, point_lat(i) - res_lat/2;
                         point_lon(i) - res_lon/2, point_lat(i) + res_lat/2;
                         point_lon(i) + res_lon/2, point_lat(i) + res_lat/2;
                         point_lon(i) + res_lon/2, point_lat(i) - res_lat/2;
                         point_lon(i) - res_lon/2, point_lat(i) - res_lat/2;
                        ];

    end

end