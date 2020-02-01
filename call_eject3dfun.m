%  This script, written by Larry G. Mastin (USGS, lgmastin@usgs.gov), reads a series of block initial
%  velocities, sizes and ejection angles from the input file:

%  input\eject3d_inputs.txt

%  It then calls the matlab function "eject3dfun", which
%  calculates block trajectories and returns the final time of flight and landing location of each block.
%  THis script then writes out a table of landing-point locations and
%  travel times to:

%  output\eject3d_output.txt

%  The inputs  that define the range and distribution of block size, initial velocity, and
%  trajectory angles are specified in Block 1 of this script

%  For a specific volcano, some other parameters may have to be set in the module eject3dfun: for example,
    % elev         = elevation of vent, in meters above sea level
    % lapse        = thermal lapse rate in the atmosphere, in degrees C per kilometer
    % wind         = the wind speed in meters per second
    % winddir_deg  = the direction toward which the wind is blowing, in degrees east of north
    % TzeroC       = the temperature at sea level, in Celsius
    % rhor         = block density, kg/m3
    % dragred      = extent of zone of reduced drag, in meters around the vent.
    % write_outfile = 'yes' if a trajectory output file is to be written for every block ejection.
    % filename_stem = stem of file names to be written out

%% BLOCK 1.  READ INPUT FILE

input_file='input\\eject3d_inputs.txt';

fid=fopen(input_file,'r');

inputlines = char(15,80);
for i=1:16
    inputline  = fgetl(fid);
    inputlines(i,1:length(inputline))=inputline;
end

% READ VARIABLES THAT CONTROL DISTRIBUTION
n_ejects=str2double(inputlines(4,26:30));
vent_lon      = str2double(inputlines(3,36:45));
vent_lat      = str2double(inputlines(3,58:65));
thetadeg_mean = str2double(inputlines(7,54:58));      %Vertical angle 
thetadeg_std  = str2double(inputlines(7,65:69));      %std
vi_mean =  str2double(inputlines(8,54:58));           %mean initial velocity, m/s
vi_std  =  str2double(inputlines(7,65:69));           %std
logdiam_mean =  str2double(inputlines(9,54:58));      %mean log block diameter, m
logdiam_std  =  str2double(inputlines(9,65:69));       %std
phimin_degrees  =  str2double(inputlines(10,54:58));   %mean log block diameter, m
phimax_degrees  =  str2double(inputlines(10,65:69));   %std
dragtypes = inputlines(12,18:57);

%Declare variables for individual trajectories
vi = NaN(n_ejects,1);
diam = NaN(n_ejects,1);
thetadeg = NaN(n_ejects,1);
phideg   = NaN(n_ejects,1);
dragtype = cell(n_ejects,1);

%Read inputs for individual trajectories
for i=1:n_ejects
    inputline = fgetl(fid);
    vi(i)       = str2double(inputline(6:12));    %vi
    diam(i)     = str2double(inputline(14:20));   %diam
    thetadeg(i) = str2double(inputline(25:28));    %thetadeg
    phideg(i)   = str2double(inputline(32:36));   %phideg
    dragtype{i} = strtrim(inputline(37:46));               %dragtype
end
fclose(fid);
%% BLOCK 4.  DECLARE OUTPUT VARIABLES AND CALL EJECT3DFUN

%Declare output variables
xfinal = zeros(n_ejects,1);
yfinal = zeros(n_ejects,1);
zfinal = zeros(n_ejects,1);
tfinal = zeros(n_ejects,1);

%Open output file, Write table header
fid = fopen('output\eject3d_output.txt','w');                %open output file

%Summarize input distributions
now = datestr(datetime);
fprintf(fid,'Inputs for eject3d, written %s\n\n',now);
fprintf(fid,'         vent location:  longitude=%10.5f,  latitude=%7.5f\n',vent_lon,vent_lat);
fprintf(fid,'      number of blocks:  %5d\n\n',n_ejects);
fprintf(fid,'                  VARIABLE            DISTRIBUTION MEAN OR MIN  STD OR MAX\n');
fprintf(fid,' vertical angle, deg. from horizontal   Gaussian      %4.1f       %4.1f\n',thetadeg_mean,thetadeg_std);
fprintf(fid,'     initial velocity, m/s              Gaussian    %6.1f      %5.1f\n',vi_mean,vi_std);
fprintf(fid,'  log10 block diameter, m               Gaussian      %4.1f       %4.1f\n',logdiam_mean,logdiam_std);
fprintf(fid,' initial direction, deg. E of N         Uniform      %5.1f      %5.1f\n\n',phimin_degrees,phimax_degrees);
fprintf(fid,'     block types:  %s\n\n',dragtypes);

%Write output table header
fprintf('   i      vi    diam   theta     phi dragtype    xfinal    yfinal    zfinal  distance    tfinal\n');
fprintf('         m/s       m     deg     deg               m          m         m        m         s\n');
fprintf(fid,'   i      vi    diam   theta     phi dragtype    xfinal    yfinal    zfinal  distance    tfinal\n');
fprintf(fid,'         m/s       m     deg     deg               m          m         m        m         s\n');

%Call eject3dfunc
for i=1:n_ejects
    dragtypenow = dragtype{i};                                   %set dragtypenow
    [xfinal(i), yfinal(i), zfinal(i), tfinal(i)] = ...
            eject3dfun(i,diam(i),vi(i),thetadeg(i),phideg(i),dragtypenow,...
            vent_lon,vent_lat);
    distance = sqrt(xfinal(i)^2 + yfinal(i)^2);
    fprintf('%4d%8.1f%8.2f%8.1f%8.1f%8s%10.1f%10.1f%10.1f%10.1f%10.1f\n', ...
            i,vi(i),diam(i),thetadeg(i),phideg(i),dragtypenow, ...
            xfinal(i),yfinal(i),zfinal(i),distance,tfinal(i));
    fprintf(fid,'%4d%8.1f%8.2f%8.1f%8.1f%8s%10.1f%10.1f%10.1f%10.1f%10.1f\n', ...
            i,vi(i),diam(i),thetadeg(i),phideg(i),dragtypenow, ...
            xfinal(i),yfinal(i),zfinal(i),distance,tfinal(i));
end

fclose(fid);
fprintf('All done\n')
