% Script that makes kml file of block trajectories

%  Written by Larry Mastin, January 2020.

% This script reads from a series of output files:
%    output\Halemaumau_run001.txt
%    output\Halemaumau_run002.txt
%    output\Halemaumau_run003.txt
%    output\Halemaumau_run004.txt
%    . . .

% . . . and writes the block trajectories to the file:
%    output\block_trajectories.kml

%  The input files give the vent location, and the x, y, and z coordinates
%  along the block trajectory.  x, y, and z are given relative to the vent.
%  This script also loads the file terrain.mat, which contains the terrain
%  coordinates and geographic projection system, to convert the x,y, and z
%  coordinates to geographic latitude, longitude and elevation asl.

%% BLOCK 1:  IMPORT MAP INFORMATION

load('input\terrain.mat');

kml_file = 'output\\block_trajectories.kml';

%% OPEN TRAJECTORY FILES AND READ THEM

n_ejects=50;
for i=1:n_ejects
    outfile = sprintf('output\\Halemaumau_run%03d.txt',i);
    fid = fopen(outfile,'r');
    if fid ~= -1
        headerlines = char(22,80);        
        % Read headerlines
        for j=1:24
            inputline = fgetl(fid);
            headerlines(j,1:length(inputline)) = inputline;
        end
        diam = str2double(headerlines(4,43:48));          %block diameter, m
        density = str2double(headerlines(5,43:48));       %block density,kg/m3
        shape = strtrim(headerlines(6,43:48));            %block shape
        vi    = str2double(headerlines(7,43:48));         %initial velocity, m/s
        thetadeg = str2double(headerlines(8,43:48));      %angle from horizontal, deg.
        phideg   = str2double(headerlines(9,43:48));      %degrees east of north
        wind     = str2double(headerlines(10,43:48));     %windspeed, m/s
        winddir  = str2double(headerlines(11,43:48));     %wind direction, deg. E of N
        Temp0z   = str2double(headerlines(16,48:52));     %temp at sea level, C
        dTdz     = str2double(headerlines(17,48:52));     %dT/dz, deg. C/km
        elev     = str2double(headerlines(18,48:52));     %vent elevation, m asl
        vent_lon = str2double(headerlines(19,49:58));     %vent longitude
        vent_lat = str2double(headerlines(20,48:56));     %vent latitude
        
        % Read trajectory
        inputline = fgetl(fid);
        clear time x y z traj_lat traj_lon;
        clear traj_x traj_y traj_z jmax;
        j=0;           %line counter in trajectory
        while ~strcmp(inputline(1:5),'#####')
            j=j+1;
            time(j) = str2double(inputline(7:13));
            x(j)    = str2double(inputline(15:23));
            y(j)    = str2double(inputline(25:33));
            z(j)    = str2double(inputline(34:43));
            inputline = fgetl(fid);
        end
   jmax = j;
   fclose(fid);
   
   %find vent UTM coordinates, vent elevation
   [vent_x, vent_y] = mfwdtran(utmstruct,vent_lat,vent_lon);           
   vent_z           = elev;
   
   %find lat, lon of each point in the trajectory
   [traj_lat, traj_lon] = minvtran(utmstruct,(vent_x+x),(vent_y+y));
   traj_z = z + elev;
   
   if i==1
       fid2 = write_kml_header(kml_file,vent_lon,vent_lat);
   end
   write_trajectory(kml_file,i,traj_lon,traj_lat,traj_z);

   end
        
end
close_kml(kml_file);   

function fid = write_kml_header(filename,vent_lon,vent_lat)

% Function that writes a kml header

fid = fopen(filename,'w');

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

% Mark vent location
fprintf(fid,'  <Placemark>\n');
fprintf(fid,'    <name>Vent Location</name>\n');
fprintf(fid,'    <styleUrl>#VentMarker</styleUrl>\n');
fprintf(fid,'    <description>vent location</description>\n');
fprintf(fid,'    <Point>\n');
fprintf(fid,'      <coordinates>%12.6e,%12.6e,0</coordinates>\n',vent_lon,vent_lat);
fprintf(fid,'    </Point>\n');
fprintf(fid,'  </Placemark>\n');

% Creat folder for trajectories
fprintf(fid,'   <Folder id="Block trajectories">\n');
fprintf(fid,'       <name>trajectories</name>\n');
fprintf(fid,'       <open>1</open>\n');
fprintf(fid,'       <visibility>1</visibility>\n');
fprintf(fid,'       <Style>\n');
fprintf(fid,'           <ListStyle>\n');
fprintf(fid,'               <bgColor>00ffffff</bgColor>\n');
fprintf(fid,'           </ListStyle>\n');
fprintf(fid,'       </Style>\n');
fclose(fid);
return
end

function write_trajectory(kml_file,j,traj_lon,traj_lat,traj_z)
% Create a trajectory
fid=fopen(kml_file,'a');
fprintf(fid,'       <Placemark>\n');
fprintf(fid,'          <name>block %d</name>\n',j);
fprintf(fid,'          <LineString>\n');
fprintf(fid,'            <extrude>0</extrude>\n');
fprintf(fid,'            <tessellate>0</tessellate>\n');
fprintf(fid,'            <altitudeMode>absolute</altitudeMode>\n');
fprintf(fid,'            <coordinates>\n');
for i=1:length(traj_lon)
    fprintf(fid,'              %12.6e,%12.6e,%6.2f\n',traj_lon(i),traj_lat(i),traj_z(i));
end
fprintf(fid,'            </coordinates>\n');
fprintf(fid,'          </LineString>\n');
fprintf(fid,'        </Placemark>\n');
fclose(fid);
end

function close_kml(kml_file)
%Close block landing points folder
fid = fopen(kml_file,'a');
fprintf(fid,'   </Folder>\n');
fprintf(fid,'</Document>\n');
fprintf(fid,'</kml>\n');
fclose(fid);

end