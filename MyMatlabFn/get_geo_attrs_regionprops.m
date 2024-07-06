%%
function E = get_geo_attrs_regionprops(lon, lat, E)

    %/ NOTE:
    %/ This function is designed for deriving the fitted ellipsed on the geographic map based on
    %/ the matlab bulit-in function regionprops().
    %/ The input E comes from E = regionprops().

    %/ *** Some original codes from Akira Agata ***
    %/ https://www.mathworks.com/matlabcentral/answers/495720-how-to-fit-an-ellipse-to-an-image-in-matlab
    
    Centr      = E.Centroid;         %/ lon, lat
    Orient     = E.Orientation;
    MajLen     = E.MajorAxisLength;
    MinLen     = E.MinorAxisLength;
    Eccen      = E.Eccentricity;
    
    %/ In this case, BW_2D has been transposed. When drawing the ellipse on geomap,
    %/ the orientation angle is in fact the angle betn minor axis and
    %/ vertical axis. Multiply by -1 to correct it.
    Orient = Orient * -1;    
    phi    = Orient;

    %/The properties are based on lon lat if given, otherwise based on the pixel size (NOTE: it starts from 1 !!!).
    dlon = abs(lon(2) - lon(1));
    dlat = abs(lat(2) - lat(1));

    %/ NOTE: regionprops define centroids by pixels. We need to get
    %        the true lon lat of the centroid.
    %/ Check https://www.mathworks.com/help/images/ref/regionprops.html
    Centr(1) = lon(floor(Centr(1))) + dlon*mod(Centr(1), 1);   
    Centr(2) = lat(floor(Centr(2))) + dlat*mod(Centr(2), 1);   

    MajLen = sqrt((MajLen.*cosd(Orient) * dlon).^2 + (MajLen.*sind(Orient) * dlat).^2 );
    MinLen = sqrt((MinLen.*sind(Orient) * dlon).^2 + (MinLen.*cosd(Orient) * dlat).^2 );

%     %/ Derive the boundary data (and other attributes) from regionprops.
%     [bndry_data_bc, MajLen_geo_bc, xMinor_bc, yMinor_bc, xMajor_bc, yMajor_bc] = ...
%             get_boundary_from_regionprops('Centr', Centr(n,:), 'MajLen', MajLen(n), 'MinLen', MinLen(n), 'phi', phi);

    theta         = linspace(0,2*pi);
    col           = (MajLen/2)*cos(theta);
    row           = (MinLen/2)*sin(theta);
    M             = makehgtform('translate',[Centr, 0],'zrotate',deg2rad(phi));
    D             = M*[col;row;zeros(1,numel(row));ones(1,numel(row))];    
    bndry_data    = D(1:2,:)';

    %/ for plotting and estimating the true geodistance.
    xMinor = Centr(1) + [ -1,  1 ] * MinLen/2*sind(phi);  %/ X-axis of two vertices of Minor axis
    yMinor = Centr(2) + [  1, -1 ] * MinLen/2*cosd(phi);  %/ Y-axis of two vertices of Minor axis
    
    xMajor = Centr(1) + [ -1, 1 ] * MajLen/2*cosd(phi);  %/ X-axis of two vertices of Major axis
    yMajor = Centr(2) + [ -1, 1 ] * MajLen/2*sind(phi);  %/ Y-axis of two vertices of Major axis

    %/ Measure the true length scale from the fitted ellipse.
    %/ distance(lat1,lon1,lat2,lon2)
    a = 6371e3; %/ Earth's radius (m)
    [ARCLEN, ~]   = distance('gc',[yMajor(1), xMajor(1)],[yMajor(2), xMajor(2)]); %/use great circle here since we are looking for shortest distance on a sphere, no matter the true course direction keeps changing.

%         one_deg_lat = a*1*pi/180*1e-3; %/ ~111 km;
    MajLen_geo = ARCLEN*pi/180*a*1e-3; %/ arc degree -> km 

    %/ Estimate the area under the ellipse.
    grid_area_2D = calc_grid_area('lon', lon, 'lat', lat)*1e-6; %/ m^2 -> km^2
    Area         = masked_by_bndry('X', grid_area_2D, 'lon', lon, 'lat', lat, 'bndry_data', bndry_data);
    
    %
    E.Centroid = Centr;         %/ lon, lat
    E.Orientation = Orient;
    E.MajorAxisLength = MajLen;
    E.MinorAxisLength = MinLen;
    E.Eccentricity = Eccen;
    E.bndry_data = bndry_data;
    E.MajLen_geo = MajLen_geo;
    E.xMinor = xMinor;
    E.yMinor = yMinor;
    E.xMajor = xMajor;
    E.yMajor = yMajor;
    E.Area = nansum(Area,'all');         %/ lon, lat
end