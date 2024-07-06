function write_nc(varargin)
    
    %/ create a set of valid parameters and their default value
    pnames = {  'ncfilename',        'data', 'data_shortname', 'data_standardname', 'data_longname', 'data_units',...
              'DeflateLevel',         'lon',            'lat',              'plev',          'time',  'time_unit',...
                      'date', 'date_format',         'remark',          'othervars'};
    
    dflts  = cell(length(pnames), 1);
    
    %/ parse function arguments
    [             ncfilename,          data,   data_shortname,   data_standardname,   data_longname,   data_units,...
                DeflateLevel,           lon,              lat,                plev,            time,    time_unit,...
                        date,   date_format,           remark,            othervars] ...
           = internal.stats.parseArgs(pnames, dflts, varargin{:});
    
    %%
    %===============================================================================
    %/      Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last Update: 15 Dec 2023
    %/
    %/ Description: This function is designed to write data into a NetCDF file
    %===============================================================================
    
    if isfile(ncfilename)
        warning('Same filename is found. Deleting it and overwriting it!') %/ This avoids the bug in overwriting a nc file.
        delete(ncfilename); 
    end
    
    if ~isempty(time)
        ntime = length(time);
    elseif ~isempty(date)
        ntime = length(date);
    else
        ntime = [];
    end
    
    nc_format = 'netcdf4_classic';  %/ A common nc format

    if isempty(DeflateLevel)    
        %/ Level of compression (integer, 0 to 9), A value of 0 indicates no compression.
        %/ A value of 1 indicates the least compression, and a value of 9 indicates the most.
        %/ You cannot set DeflateLevel if Datatype is "string".
        DeflateLevel = 2;   %/ by default, setting 2 is good enough.
    end  

    %/ Always convert cell into char
    if iscell(data_shortname)     
        data_shortname    = char(data_shortname);      
    end
    if iscell(data_standardname)  
        data_standardname = char(data_standardname);   
    end
    if iscell(data_longname)      
        data_longname     = char(data_longname);       
    end
    if iscell(data_units)         
        data_units        = char(data_units);          
    end
    if iscell(remark)             
        remark            = char(remark);              
    end

    %/ If main data is given
    if ~isempty(data)
        data_vtype = class(data);   %/ Check data type
        ndim = length(size(data));
    
        %/ Define data dimensions
        if ~isempty(ntime) && ~isempty(lon) && ~isempty(lat)
            if ~isempty(plev)
                Dimensions = {'lon', length(lon), 'lat', length(lat), 'plev', length(plev), 'time', ntime};
            else
                Dimensions = {'lon', length(lon), 'lat', length(lat), 'time', ntime};
            end
            
        elseif ~isempty(ntime) && isempty(lon) && isempty(lat)
            Dimensions = {'time', ntime};
            
        elseif isempty(ntime) && ~isempty(lon) && ~isempty(lat)
            Dimensions = {'lon', length(lon), 'lat', length(lat)};
        
        elseif isempty(ntime) && isempty(lon) && isempty(lat) %/ Otherwise, define dimension name arbitrarily
            if ndim == 2
                Dimensions = {'dim1', size(data, 1), 'dim2', size(data, 2)};
            elseif ndim == 3
                Dimensions = {'dim1', size(data, 1), 'dim2', size(data, 2), 'dim3', size(data, 3)};
            else
                error('code not set for ndim == %d!', ndim);
            end
        else
            error('code not set!');
        end
        nccreate(ncfilename, data_shortname, 'datatype', data_vtype, 'Dimensions', Dimensions, 'DeflateLevel', DeflateLevel, 'format', nc_format);
        
        ncwriteatt(ncfilename, data_shortname, 'standard_name', data_standardname);
        ncwriteatt(ncfilename, data_shortname, 'long_name',     data_longname);
        ncwriteatt(ncfilename, data_shortname, 'units',         data_units);
        ncwriteatt(ncfilename, data_shortname, 'remark',        remark);
        % ncwriteatt(ncfilename, data_shortname, 'missing_value', '-1e+04');
        
        if ~isempty(lon)
            lon_vtype = class(lon);
            nccreate(ncfilename,   'lon', 'datatype', lon_vtype, 'Dimensions',{'lon', length(lon)}, 'format', nc_format);
            ncwriteatt(ncfilename, 'lon', 'standard_name', 'longitude');
            ncwriteatt(ncfilename, 'lon', 'long_name', 'longitude');
            ncwriteatt(ncfilename, 'lon', 'units', 'degrees_east');
            ncwriteatt(ncfilename, 'lon', '_CoordinateAxisType', 'lon');
            ncwrite(ncfilename,'lon',lon);
        end
        
        if ~isempty(lat)
            lat_vtype = class(lat);
            nccreate(ncfilename,   'lat', 'datatype', lat_vtype, 'Dimensions',{'lat', length(lat)}, 'format', nc_format);
            ncwriteatt(ncfilename, 'lat', 'standard_name', 'latitude');
            ncwriteatt(ncfilename, 'lat', 'long_name', 'latitude');
            ncwriteatt(ncfilename, 'lat', 'units', 'degrees_north');
            ncwriteatt(ncfilename, 'lat', '_CoordinateAxisType', 'lat');
            ncwrite(ncfilename,'lat',lat);
        end
        
        if ~isempty(plev)
            plev_vtype = class(plev);
            nccreate(ncfilename,   'plev', 'datatype', plev_vtype, 'Dimensions',{'plev', length(plev)}, 'format', nc_format);
            ncwriteatt(ncfilename, 'plev', 'standard_name', 'plev');
            ncwriteatt(ncfilename, 'plev', 'long_name', 'pressure level');
            ncwriteatt(ncfilename, 'plev', 'units', 'Pa');
            ncwriteatt(ncfilename, 'plev', '_CoordinateAxisType', 'plev');
            ncwrite(ncfilename,'plev',plev);
        end
        
        if ~isempty(time)
            time_vtype = class(time);
            nccreate(ncfilename,   'time', 'datatype', time_vtype, 'Dimensions',{'time', length(time)}, 'format', nc_format);
            ncwriteatt(ncfilename, 'time', 'long_name', 'time');
            
            if isempty(time_unit)
                ncwriteatt(ncfilename, 'time', 'units', 'hours since 1900-01-01 00:00:00.0'); %/ by default
            else
                ncwriteatt(ncfilename, 'time', 'units', time_unit);
            end
            ncwriteatt(ncfilename, 'time', '_CoordinateAxisType', 'time');
            ncwrite(ncfilename,'time',time);
        end
        
        if ~isempty(date)
            date_vtype = class(date);
            if isequal(date_vtype, 'int64')
                date = int32(date);  %/ For 'int64' is not a supported datatype for a netcdf4_classic format file. Use int32 instead.
                date_vtype = 'int32';
            end
    
            nccreate(ncfilename,   'date', 'Datatype', date_vtype, 'Dimensions',{'date', length(date)}, 'format', nc_format); %/ sadly int64 is somehow not supported...
            ncwriteatt(ncfilename, 'date', 'long_name', 'date');
            
            if ~isempty(date_format)
                ncwriteatt(ncfilename, 'date', 'units', date_format);
            else
                ncwriteatt(ncfilename, 'date', 'units', 'yyyymmdd');
            end
            ncwriteatt(ncfilename, 'date', '_CoordinateAxisType', 'date');
            ncwrite(ncfilename,'date',date);
        end
        % ncdisp(ncfilename)
        ncwrite(ncfilename,data_shortname,data);
    end

    %/ Loop over the fields in the struct 'othervars' and write into nc
    if ~isempty(othervars)
        if ~isstruct(othervars)
            error('othervars must be a struct!');
        end
        
        fld = fieldnames(othervars);
        for k = 1:length(fld)
            fprintf('Storing the field %s into the nc file...\n', fld{k});
            vname  = char(fld(k));
            var    = othervars.(vname);
            if iscell(var)
                var = char(var);  %/ Convert cell to char
            end
            vtype  = class(var);
            if isequal(vtype, 'int64')
                var = int32(var);  %/ For 'int64' is not a supported datatype for a netcdf4_classic format file. Use int32 instead.
                vtype = 'int32';
            end

            vshape = size(var);
            ndims  = length(vshape);
            vlen   = length(var(:));
            switch vtype
                case {'double','single','int32'}
                    if vlen==1
                        nccreate(ncfilename, vname, 'Datatype',vtype,'format',nc_format);

                    elseif contains(vname, '_lon') || contains(vname, '_lat') 
                        if contains(vname, '_lon')
                            vlongname = 'longitude';
                            vunits    = 'degrees_east';
                        else
                            vlongname = 'latitude';
                            vunits    = 'degrees_north';
                        end
                        nccreate(ncfilename,   vname, 'datatype', vtype, 'Dimensions',{vname, length(var)}, 'format', nc_format);
                        ncwriteatt(ncfilename, vname, 'standard_name', vlongname);
                        ncwriteatt(ncfilename, vname, 'long_name',     vlongname);
                        ncwriteatt(ncfilename, vname, 'units', vunits);
                        ncwriteatt(ncfilename, vname, '_CoordinateAxisType', vname);
                        ncwrite(ncfilename,vname,var);
                        
                    else
                        vname_lon = strcat(vname, '_lon');
                        vname_lat = strcat(vname, '_lat');
                        if isfield(othervars, vname_lon) && isfield(othervars, vname_lat)
                            Dimensions_othervar = {vname_lon, length(othervars.(vname_lon)), vname_lat, length(othervars.(vname_lat))};

                        elseif isequal(size(data), size(var)) %/ If var's dims == the main data's, then assume it shares the same lon/lat/lev/time
                            Dimensions_othervar = Dimensions;

                        else
                            if min(vshape)==1
                                Dimensions_othervar = {[vname '1'] vlen};  %/ Instead of a common name dim1 and dim2 because they must be *unique*!
                            elseif ndims==2
                                Dimensions_othervar = {[vname '1'] vshape(1) [vname '2'] vshape(2)};
                            elseif ndims==3
                                Dimensions_othervar = {[vname '1'] vshape(1) [vname '2'] vshape(2) [vname '3'] vshape(3)};
                            else
                                error('code not set for ndims == %d!', ndims)
                            end
                        end
                        nccreate(ncfilename, vname, 'Datatype', vtype, 'Dimensions', Dimensions_othervar, 'DeflateLevel', DeflateLevel,  'format',nc_format);
                    end

                case {'char'}
                    if ndims == 2
                        Dimensions_othervar = {[vname '1'] vshape(1) [vname '2'] vshape(2)};
                    else
                        error('code not set for ndims == %d!', ndims)
                    end
                    nccreate(ncfilename, vname, 'Datatype', vtype, 'Dimensions',Dimensions_othervar, 'format',nc_format);

                otherwise
                    error('Detected othervars.%s is a %s! Only ''double'', ''single'' or ''char'' are valid!', vname, vtype);
            end
            ncwrite(ncfilename,vname,var); %/ Write data
        end
    end
    fprintf('!!! NetCDF data is written into %s !!!\n', ncfilename)
    ncdisp(ncfilename)
end