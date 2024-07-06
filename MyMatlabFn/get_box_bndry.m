function bndry_data = get_box_bndry(lon, lat)

    if length(lon) ~= 2   error('Input lon must a two-element array!');  end
    if length(lat) ~= 2   error('Input lat must a two-element array!');  end

    bndry_data = [lon(1), lat(1);
                  lon(1), lat(2);
                  lon(2), lat(2);
                  lon(2), lat(1);
                  lon(1), lat(1);];
end