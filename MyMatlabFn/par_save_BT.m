%%
function par_save_BT(fname, lon, lat)
    BT = [];
    BT.lon_time = lon;
    BT.lat_time = lat;
    save(fname, 'BT', '-v7.3')
  
end

