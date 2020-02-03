%Script that reads the Kilauea GEOTIFF file and writes out the file
%     terrain.mat containing:
%     terrain  %An NxM matrix of terrain elevations, meters asl
%     t_lat    %An NxM matrix of latitudes for each node in terrain
%     t_lon    %An NxM matrix of longitudes for each node in terrain
%     t_x      %An NxM matrix of UTM x coordinates for each node in terrain
%     t_y      %an NxM matrix of UTM y coordinates for each node in terrain
%     proj     %variables extracted from geotiff giving projection values'
%     geotiff_name   Name of geotiff file that was read.

%%  Variables to be chcked and/or changed before running:


% Vent longitude, latitude
vent_lon = -155.2814;                     %Lat, lon
vent_lat = 19.4080;

%Name of GEOTIFF file.  It is assumed that this file has coordinates in
%lat/lon, and that the values of each pixel are the elevation asl in
%meters.
geotiff_name = 'input\KIL_30m2019\KIL_30m2019.tif';

%% Start calculations

fprintf('Running tifreader\n');
fprintf('vent_lon  = %9.4f\n',vent_lon);
fprintf('vent_lat  = %9.4f\n',vent_lat);

%%  Set up UTM structure to calculate distances
UTM_zone = utmzone([vent_lat vent_lon]);               
[ellipsoid,estr] = utmgeoid(UTM_zone);        %find UTM zone
utmstruct = defaultm('utm');
utmstruct.zone = UTM_zone;                    %add zone to utm structure
utmstruct.geoid = ellipsoid;            %add ellipsoid to utm structure
utmstruct = defaultm(utmstruct);
[vent_x, vent_y, UTM_zone] = deg2utm(vent_lat,vent_lon);

fprintf('UTM_ZONE  = %s\n',UTM_zone);
fprintf('vent_x    = %9.0f\n',vent_x);
fprintf('vent_y    = %9.0f\n\n',vent_y);

%% Read geotiff, get coordinates
fprintf('Reading geotiff file %s\n',geotiff_name);

%--------------------------------------------------------------------------
%Note:  you will need the Matlab mapping toolbox in order for these
%functions to work.
[terrain,refmat,bbox] = geotiffread(geotiff_name);
proj = geotiffinfo(geotiff_name);
%--------------------------------------------------------------------------

%Calculate locations of cell centers (assume bbox gives the cell centers
% at the corners)
fprintf('Determining cell locations\n');

latmin = bbox(1,2);
lonmin = bbox(1,1);
latmax = bbox(2,2);
lonmax = bbox(2,1);

dlat = (latmin-latmax)/(size(terrain,1));
dlon = (lonmax-lonmin)/(size(terrain,2));

fprintf('Bounding box values:\n');
fprintf('latitude:  %9.4f to  %9.4f\n',latmin,latmax);
fprintf('longitude: %9.4f to  %9.4f\n',lonmin,lonmax);
fprintf('dlat= %10.5e,  dlon=%10.5e\n',dlat,dlon);

%Calculate cell centers.
% Bounding box values must be offset by one half a cell width
% to get cell centers
latCC = latmax+dlat/2:dlat:latmin;      %cell centers in lat
lonCC = lonmin+dlon/2:dlon:lonmax;      %in lon

%Create meshgrid of lat, lon
[t_lon, t_lat] = meshgrid(lonCC,latCC);
[t_lon2, t_lat2] = meshgrid(lonCC,flipud(latCC));
%Transform to UTM coordinates
[t_x, t_y, UTM_zone] = deg2utm(t_lat,t_lon);

%Write out projection reaults
fprintf('saving results to terrain.mat\n');
save('input\terrain.mat','terrain','t_lon','t_lat','t_x','t_y', ...
     'UTM_zone');
fprintf('All done\n');

