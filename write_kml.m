% Script that writes kml output file

%% Block 1:  Read outputs

output_file='output\eject3d_output.txt';

fid=fopen(output_file,'r');

inputlines = char(15,80);
for i=1:15
    inputline  = fgetl(fid);
    inputlines(i,1:length(inputline))=inputline;
end

% READ VARIABLES THAT CONTROL DISTRIBUTION
vent_lon      = str2double(inputlines(3,36:45));
vent_lat      = str2double(inputlines(3,58:65));
n_ejects=str2double(inputlines(4,26:30));
thetadeg_mean = str2double(inputlines(7,54:58));      %Vertical angle 
thetadeg_std  = str2double(inputlines(7,65:69));      %std
vi_mean =  str2double(inputlines(8,54:58));           %mean initial velocity, m/s
vi_std  =  str2double(inputlines(8,65:69));           %std
logdiam_mean =  str2double(inputlines(9,54:58));      %mean log block diameter, m
logdiam_std  =  str2double(inputlines(9,65:69));       %std
phimin_degrees  =  str2double(inputlines(10,54:58));   %mean log block diameter, m
phimax_degrees  =  str2double(inputlines(10,65:69));   %std
dragtype = inputlines(12,18:57);

%Declare variables for individual trajectories
vi = NaN(n_ejects,1);
diam = NaN(n_ejects,1);
thetadeg = NaN(n_ejects,1);
phideg   = NaN(n_ejects,1);
dragtype = cell(n_ejects,1);
xfinal   = NaN(n_ejects,1);
yfinal   = NaN(n_ejects,1);
zfinal   = NaN(n_ejects,1);
distance = NaN(n_ejects,1);
tfinal   = NaN(n_ejects,1);

%Read inputs for individual trajectories
for i=1:n_ejects
    inputline = fgetl(fid);
    vi(i)       = str2double(inputline(6:12));    %vi
    diam(i)     = str2double(inputline(14:20));   %diam
    thetadeg(i) = str2double(inputline(25:28));    %thetadeg
    phideg(i)   = str2double(inputline(32:36));   %phideg
    dragtype{i} = strtrim(inputline(37:46));      %dragtype
    xfinal(i)   = str2double(inputline(47:54));
    yfinal(i)   = str2double(inputline(57:64));
    zfinal(i)   = str2double(inputline(67:74));
    distance(i) = str2double(inputline(77:84));
    tfinal(i)   = str2double(inputline(87:94));
end
fclose(fid);

%% BLOCK 2:  LOAD TERRAIN AND MAP DATA, FIND GEOGRAPHIC LOCATIONS

load('input\terrain.mat','terrain','t_lon','t_lat','utmstruct');                  %load terrain file

[vent_x, vent_y]       = mfwdtran(utmstruct,vent_lat,vent_lon);                   %get vent coordinates in UTM
vent_z                 = double(interp2(t_lon,t_lat,terrain,vent_lon,vent_lat));  %get vent elevation (m)
[lat_final, lon_final] = minvtran(utmstruct,(vent_x+xfinal),(vent_y+yfinal));     %get landing point coordinates, lat/lon
elev_final             = zfinal + vent_z;                                         %get landing point elevation

%% BLOCK 3:  START WRITING KML FILE

fid = fopen('output\eject3d_output.kml','w');

fprintf(fid,'<?xml version="1.0" encoding="UTF-8"?>\n');
fprintf(fid,'<kml xmlns="http://www.opengis.net/kml/2.2">\n');
fprintf(fid,'<Document>\n');

% Vent location marker
fprintf(fid,'    <StyleMap id="VentMarker">\n');
fprintf(fid,'      <Pair>\n');
fprintf(fid,'             <key>normal</key>\n');
fprintf(fid,'             <styleUrl>#VentMarkerNormal</styleUrl>\n');
fprintf(fid,'      </Pair>\n');
fprintf(fid,'      <Pair>\n');
fprintf(fid,'             <key>highlight</key>\n');
fprintf(fid,'             <styleUrl>#VentMarkerHighlight</styleUrl>\n');
fprintf(fid,'      </Pair>\n');
fprintf(fid,'    </StyleMap>\n');
fprintf(fid,'    <Style id="VentMarkerNormal">\n');
fprintf(fid,'      <IconStyle>\n');
fprintf(fid,'             <scale>0.8</scale>\n');
fprintf(fid,'             <Icon>\n');
fprintf(fid,'               <href>http://maps.google.com/mapfiles/kml/pal4/icon52.png</href>\n');
fprintf(fid,'             </Icon>\n');
fprintf(fid,'      </IconStyle>\n');
fprintf(fid,'      <LabelStyle>\n');
fprintf(fid,'             <scale>0</scale>\n');
fprintf(fid,'      </LabelStyle>\n');
fprintf(fid,'    </Style>\n');
fprintf(fid,'    <Style id="VentMarkerHighlight">\n');
fprintf(fid,'      <IconStyle>\n');
fprintf(fid,'             <scale>1</scale>\n');
fprintf(fid,'             <Icon>\n');
fprintf(fid,'               <href>http://maps.google.com/mapfiles/kml/pal4/icon52.png</href>\n');
fprintf(fid,'             </Icon>\n');
fprintf(fid,'      </IconStyle>\n');
fprintf(fid,'      <LabelStyle>\n');
fprintf(fid,'             <scale>1</scale>\n');
fprintf(fid,'      </LabelStyle>\n');
fprintf(fid,'    </Style>\n');

%Block markers
fprintf(fid,'    <StyleMap id="BlockMarker">\n');
fprintf(fid,'      <Pair>\n');
fprintf(fid,'             <key>normal</key>\n');
fprintf(fid,'             <styleUrl>#BlockMarkerNormal</styleUrl>\n');
fprintf(fid,'      </Pair>\n');
fprintf(fid,'      <Pair>\n');
fprintf(fid,'             <key>highlight</key>\n');
fprintf(fid,'             <styleUrl>#BlockMarkerHighlight</styleUrl>\n');
fprintf(fid,'      </Pair>\n');
fprintf(fid,'    </StyleMap>\n');
fprintf(fid,'    <Style id="BlockMarkerNormal">\n');
fprintf(fid,'      <IconStyle>\n');
fprintf(fid,'             <scale>0.8</scale>\n');
fprintf(fid,'             <Icon>\n');
fprintf(fid,'               <href>http://maps.google.com/mapfiles/kml/pal4/icon48.png</href>\n');
fprintf(fid,'             </Icon>\n');
fprintf(fid,'      </IconStyle>\n');
fprintf(fid,'      <LabelStyle>\n');
fprintf(fid,'             <scale>0</scale>\n');
fprintf(fid,'      </LabelStyle>\n');
fprintf(fid,'    </Style>\n');
fprintf(fid,'    <Style id="BlockMarkerHighlight">\n');
fprintf(fid,'      <IconStyle>\n');
fprintf(fid,'             <scale>1</scale>\n');
fprintf(fid,'             <Icon>\n');
fprintf(fid,'               <href>http://maps.google.com/mapfiles/kml/pal4/icon48.png</href>\n');
fprintf(fid,'             </Icon>\n');
fprintf(fid,'      </IconStyle>\n');
fprintf(fid,'      <LabelStyle>\n');
fprintf(fid,'             <scale>1</scale>\n');
fprintf(fid,'      </LabelStyle>\n');
fprintf(fid,'    </Style>\n');

% Mark vent location
fprintf(fid,'  <Placemark>\n');
fprintf(fid,'    <name>Vent Location</name>\n');
fprintf(fid,'    <styleUrl>#VentMarker</styleUrl>\n');
fprintf(fid,'    <description><![CDATA[\n');
fprintf(fid,'                <ul>\n');
fprintf(fid,'                  <li>Number of blocks ejected: %d</li>\n',n_ejects);
fprintf(fid,'                  <li>     vertical angle: mean=%4.1f, std=%4.1d deg from horizontal</li>\n',...
    thetadeg_mean,thetadeg_std);
fprintf(fid,'                  <li>   horizontal angle: from %6.1f to %6.1f deg E of N</li>\n',...
    phimin_degrees,phimax_degrees);
fprintf(fid,'                  <li> log block diameter: mean=%4.2f, std=%4.2f meters</li>\n',...
    logdiam_mean,logdiam_std);
fprintf(fid,'                </ul>]]>\n');
fprintf(fid,'    </description>\n');
fprintf(fid,'    <Point>\n');
fprintf(fid,'      <coordinates>%12.6e,%12.6e,0</coordinates>\n',vent_lon,vent_lat);
fprintf(fid,'    </Point>\n');
fprintf(fid,'  </Placemark>\n');

% Create folder of blocks
fprintf(fid,'   <Folder id="block landing points">\n');
fprintf(fid,'       <name>Forecasts</name>\n');
fprintf(fid,'       <open>1</open>\n');
fprintf(fid,'       <visibility>1</visibility>\n');
fprintf(fid,'       <Style>\n');
fprintf(fid,'           <ListStyle>\n');
fprintf(fid,'               <bgColor>00ffffff</bgColor>\n');
fprintf(fid,'           </ListStyle>\n');
fprintf(fid,'       </Style>\n');

% Create landing point markers
for i=1:n_ejects
    fprintf(fid,'  <Placemark>\n');
    fprintf(fid,'    <name>Block %d</name>\n',i);
    fprintf(fid,'    <styleUrl>#BlockMarker</styleUrl>\n');
    fprintf(fid,'    <description><![CDATA[\n');
    fprintf(fid,'                <ul>\n');
    fprintf(fid,'                  <li>Ejection velocity: %6.1f m/s</li>\n',vi(i));
    fprintf(fid,'                  <li>   vertical angle: %4.1f deg from horizontal</li>\n',thetadeg(i));
    fprintf(fid,'                  <li> horizontal angle: %6.1f deg E of N</li>\n',phideg(i));
    fprintf(fid,'                  <li>   block diameter: %4.1f meters</li>\n',diam(i));
    fprintf(fid,'                  <li>      block shape: %s</li>\n',dragtype{i});
    fprintf(fid,'                  <li>   final distance: %7.1f meters</li>\n',distance(i));
    fprintf(fid,'                  <li>   time in flight: %6.1f seconds</li>\n',tfinal(i));
    fprintf(fid,'                </ul>]]>\n');
    fprintf(fid,'    </description>\n');
    fprintf(fid,'    <Point>\n');
    fprintf(fid,'      <coordinates>%12.6e,%12.6e,0</coordinates>\n',lon_final(i),lat_final(i));
    fprintf(fid,'    </Point>\n');
    fprintf(fid,'  </Placemark>\n');
end


%Close block landing points folder
fprintf(fid,'   </Folder>\n');
fprintf(fid,'</Document>\n');
fprintf(fid,'</kml>\n');

fclose(fid);
fprintf('all done\n');

