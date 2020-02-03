
function [xfinal,yfinal,zfinal,tfinal] = eject3dfun(inow,diam,vi,thetadeg,phideg,dragtype,...
          vent_lon,vent_lat)

    % This function runs ballistic a trajectory calculation.  It it called
    % by the script call_eject3dfun.m.
    % inputs:
    % inow      = run number (number in a series of ejections run by
    %               call_eject3dfun.m)
    % diam      = block diameter, meters
    % thetadeg  = vertical trajectory angle, degrees from horizontal
    % phideg    = horizontal trajectory angle, degrees CW from north
    % dragttype = 'lo_cube','hi_cube','sphere','drag_shell'
    % vent_lon  = vent longitude
    % vent_lat  = vent latitude
    %
    % In addition to the call arguments, eject3dfun.m includes atmospheric
    % properties and other inputs specified in Block 1 of the code.
    %
    % If, in block 1, write_outfile='yes', then a trajectory output table
    % is written to the file:
    %   output\${filename_stem}_###.txt,
    %      where %{filename_stem} is a partiial file name given in block 1
    %      ### is the run number (inow)

    
    %% BLOCK 1.  VOLCANO-SPECIFIC PARAMETERS.
    % These parameters should be reviewed for each series of runs

    dragred = 0;                     %zone of reduced drag, meters from vent
    %elev = 1000;                     %elevation of ejection point
    lapse = 6.5;                     %thermal lapse rate, deg. K/km
    rhor = 2500;                     %density of ballistic, kg/m3
    wind = 0.;                       %wind velocity, m/s
    winddir_deg = 0.;                %direction toward which wind is blowing, degrees CW from north
    xi = 0;                          %vertical distance of landing point above ejection point, meters
                                     %  Note: Only negative values (<0) are allowed for xi.
    TzeroC = 25.;                     %Temperature at sealevel, Celsius
    write_outfile = 'yes';           %='yes' if an output file is to be written for each run:
                                    %if write_outfile='yes', the output files have the format:
                                    % eject_out_YY-MM-DD_HH_mm_run%% %.txt, where YY, MM, DD, HH and mm
                                    % are the year, month, day, hour, and minute of the run series
    filename_stem = 'output/Halemaumau';    %If write_outfile='yes', this gives the stem of filename to be 
                                    %written out.  The full filename will be:
                                    % "filename_stem_run%% %.txt", where %% % is the run number (inow)
    terrain_file = 'input\\terrain.mat';   %name of file containing digital terrain                                    

    %% BLOCK 2.  PHYSICAL CONSTANTS 

    pi   = 3.14159;                 %pi
    grav = 9.80665;                 %gravitational constant, m/s2
    R_air =286.98;                  %Specific gas constant for air (J/kg K)

    %% BLOCK 3.  PERFORM INITIAL CALCULATIONS AND ERROR CHECKS
    
    load(terrain_file,'terrain','t_lon','t_lat','utmstruct');            %load terrain file
    [vent_x, vent_y, UTM_zone] = deg2utm(vent_lat,vent_lon);
    %[vent_x, vent_y] = mfwdtran(utmstruct,vent_lat,vent_lon);            %vent location in UTM
    vent_z  = double(interp2(t_lon,t_lat,terrain,vent_lon,vent_lat));  %get vent elevation
    elev = vent_z;
    xi = -100;                   %subtract a little to make sure ejection starts
    
    %dragtype = 'shell'             %can choose 'hicube', 'locube', 'sphere', 'shell', or 'const'
    if strcmp(dragtype,'const')
        cd_const   = 1;              %if dragtype='const', specify value of Cd
    end

    rad     = diam / 2;                             %convert block diameter to radius, m
    theta   = pi*thetadeg/180;                      %convert theta to radians
    phi2    = pi*(0.5-phideg/180);                  %convert phi to radians, degrees CCW from east
    xwind   = wind*cos(pi*(0.5-winddir_deg/180));   %x (east) component of wind
    ywind   = wind*sin(pi*(0.5-winddir_deg/180));   %y (north) component of wind

    TzeroK  = TzeroC + 273.15;                      %convert temperature to Kelvin
%    ttotcd0 = (-vi * sin(theta) - ...
%               sqrt(vi ^ 2 * (sin(theta) ^ 2) - ...
%               2 * grav * xi))/ (-grav);            %estimate time of travel, assuming cd=0
    ttotcd0 = (-vi * sin(theta) - ...
               sqrt(vi ^ 2 * (sin(theta) ^ 2)))/ (-grav);            %estimate time of travel, assuming cd=0
    dt      = ttotcd0/100;                          %Give initial time step as a small fraction of total time

    %%  BLOCK 4.  DECLARE OUTPUT ARRAYS

    time=NaN(5000,1); 
    x=NaN(5000,1);                  %vector for horizontal position from vent (m)
    y=NaN(5000,1);                  %vector for horizontal position from vent (m)
    z=NaN(5000,1);                  %vector for vertical position above vent (m)
    vx=NaN(5000,1);                %vector for x-(east) components of velocity (m/s)
    vy=NaN(5000,1);                %vector for y-(north) components of velocity (m/s)
    vz=NaN(5000,1);                %vector for z-components of velocity (m/s)
    mach=NaN(5000,1);            %mach number now
    temp=NaN(5000,1);            %air temperature now
    rhoa=NaN(5000,1);            %air density now
    reynolds = NaN(5000,1);  %Reynolds number now
    cdnow=NaN(5000,1);          %cd now

    %Provide values to first elements in arrays
    x(1)    = 0.;                                   %initialize x
    y(1)    = 0.;                                   %initialize y
    z(1)    = 0.;                                   %initialize z
    vx(1)   = vi*cos(theta)*cos(phi2);              %initial vx
    vy(1)   = vi*cos(theta)*sin(phi2);              %initial vy
    vz(1)   = vi*sin(theta);                        %initial vz
    %v=vi                                           %initial velocity, m/s
    time(1) = 0;                                    %initial time (seconds)
    %zmax    = 0 ;                                   %maximum point in block's trajectory
    po = 101300 * ((TzeroK - z(1) * lapse / 1000)/ TzeroK) ^ ...
        (-grav / (R_air ^ 2 * lapse / 1000));  %air pressure at vent

    %%  BLOCK 5.  CALCULATE BLOCK MASS AND FRONTAL AREA DEPENDING ON DRAG TYPES

    if strcmp(dragtype,'hicube')
        mass = 8 * rad ^ 3 * rhor;  %mass of clast calculated based on the option "high cube"
        area = 4 * rad * rad;       %area of clast calculated based on the option "high cube"
    elseif strcmp(dragtype,'locube')
        mass = 8 * rad ^3 * rhor;   %mass of a cube
        area = 4 * rad ^2 * sqrt(3);   %frontal area of a cube pointing forward
    elseif strcmp(dragtype,'sphere')
        mass = (4/3)*pi*rad^3*rhor;   %mass of a sphere
        area = pi*rad^2;              %frontal area of a sphere
    elseif strcmp(dragtype,'shell')
        mass = pi * 2.4107 * rad^3 * rhor;   %mass of a G1 shaped ogive cylinder
        area = pi * rad^2;                   %frontal area of an ogive cylinder
    elseif strcmp(dragtype,'const')
        mass = (4/3)*pi*rad^3*rhor;   %assume the mass of a sphere
        area = pi*rad^2;              %and the frontal area of a sphere
    end

    %%  BLOCK 6.  START INTEGRATING THROUGH TIME

    xnow=x(1);
    ynow=y(1);
    znow=z(1);
    vxnow=vx(1);
    vynow=vy(1);
    vznow=vz(1);
    i=1;

    %fprintf('%4d    %4.2f    %6.1f    %6.1f    %6.1f  %10.4f  %10.4f\n',i,time(i),distance,ground_z_now,(znow+elev),lat_now,lon_now))
    %fprintf('   i    time  distance      elev      znow         lat         lon\n');
    while znow >= xi %run while block elevation is above the ground elevation (xi) 
    
        %Specify current properties
        x(i)=xnow; 
        y(i)=ynow; 
        z(i)=znow;
        vx(i)=vxnow;
        vy(i)=vynow;
        vz(i)=vznow;
        temp(i) = TzeroK - (znow + elev) * lapse / 1000;              %temperature at znow
        visc = (0.000172 * (390 / (temp(i) + 117)) * ...
            (temp(i) / 273) ^ 1.5) / 10;                         %air viscosity at znow
        pressure = po *((TzeroK - (znow + elev) * ...
            lapse / 1000) / TzeroK) ^ ...
            (grav / (R_air * lapse / 1000));                      %air pressure at znow
        rhoa(i) = pressure / (R_air * temp(i));                       %air density at znow
        c_sound    = 20.116 * sqrt(temp(i));                          %sound speed at znow
        vwind   = sqrt((vxnow-xwind)^2 + (vynow-ywind)^2 ...
                + vznow^2);                                      %velocity minus wind component
        reynolds(i) = rhoa(i) * vwind * 2 * rad / visc;               %Reynolds number
        mach(i) = vwind / c_sound;                                    %Mach number

        %Calculate drag coefficient
        if strcmp(dragtype,'hicube')
            cd = drag_hicube(mach(i));
        elseif strcmp(dragtype,'locube')
            cd = drag_locube(mach(i));
        elseif strcmp(dragtype,'sphere')
            cd = drag_sphere(reynolds(i),mach(i));
        elseif strcmp(dragtype,'shell')
             cd = drag_shell(mach(i));
        elseif strcmp(dragtype,'const')
            cd=cd_const;
        end

        %calculate drag reduction in reduced-drag zone
        if (dragred>0) && (sqrt(xnow^2 + znow^2) < dragred)
            cdnow(i) = cd * (sqrt(xnow^2 + znow^2)/dragred)^2;
        else
            cdnow(i)=cd;
        end

        %%  Runge Kutta to determine new position and velocities
           
        RK=rk3d(xnow,ynow,znow,vxnow,vynow,vznow,xwind,ywind,dt,rhoa(i),cdnow(i),area,mass,grav,rhor);

        %update position and velocity
        xnow  = RK(1);
        ynow  = RK(2);
        znow  = RK(3);
        vxnow = RK(4);
        vynow = RK(5);
        vznow = RK(6);
                                       
        %update i so that next value is stored in the next position of storage vector
        i=i+1;

        %update time
        time(i) = time(i-1) + dt;
        
        %Determine current longitude, latitude, elevation
        [lat_now, lon_now]     = utm2deg((vent_x+xnow),(vent_y+ynow),UTM_zone);
        %[lat_now, lon_now]    = minvtran(utmstruct,(vent_x+xnow),(vent_y+ynow));
        ground_z_now = double(interp2(t_lon,t_lat,terrain,lon_now,lat_now));
        xi = ground_z_now-elev;
        
        distance = sqrt(xnow^2+ynow^2);

        %fprintf('%4d    %4.2f    %6.1f    %6.1f    %6.1f  %10.4f  %10.4f\n',i,time(i),distance,ground_z_now,(znow+elev),lat_now,lon_now);
        %fprintf('distance=%f, ground_z_now=%f, xi=%f, znow=%f\n',distance,ground_z_now,xi,znow);
                        
       %END OF ITERATIONS

    end

    %% BLOCK 7:  INTERPOLATE TO FIND THE FINAL VALUES

    xfinal = x(i-1) + (xnow - x(i-1)) * (xi - z(i-1)) / (znow - z(i-1));           %final x
    yfinal = y(i-1) + (ynow - y(i-1)) * (xi - z(i-1)) / (znow - z(i-1));           %final y
    zfinal = xi;                                                                   %final z
    tfinal = time(i-1) + dt * (xi - z(i-1)) / (znow - z(i-1));                          %final time (s)
    vxfinal = vx(i-1) + (vxnow - vx(i-1)) * (xi - z(i-1)) / (znow - z(i-1));       %final vx
    vyfinal = vy(i-1) + (vynow - vy(i-1)) * (xi - z(i-1)) / (znow - z(i-1));       %final vy
    vzfinal = vz(i-1) + (vznow - vz(i-1)) * (xi - z(i-1)) / (znow - z(i-1));       %final vz
    %vfinal = sqrt(vxfinal ^ 2 + vzfinal ^ 2);                                    %final velocity

    %Make these the final values of the arrays
    x(i)    = xfinal; 
    y(i)    = yfinal; 
    z(i)    = zfinal; 
    vx(i)   = vxfinal; 
    vy(i)   = vyfinal; 
    vz(i)   = vzfinal; 
    time(i) = tfinal;

    %Calculate other properties at the final elevation
    distance = sqrt(xfinal^2 + yfinal^2);
    temp(i) = TzeroK - (zfinal + elev) * lapse / 1000;              %temperature at zfinal
    visc = (0.000172 * (390 / (temp(i) + 117)) * ...
            (temp(i) / 273) ^ 1.5) / 10;                     %air viscosity at zfinal
    pressure = po *((TzeroK - (zfinal + elev) * ...
         lapse / 1000) / TzeroK) ^ ...
         (grav / (R_air * lapse / 1000));                      %air pressure at zfinal
    rhoa(i) = pressure / (R_air * temp(i));                        %air density at zfinal
    c_sound    = 20.116 * sqrt(temp(i));                          %sound speed at zfinal
    %velocity = sqrt(vxnow*vxnow + vynow*vynow + vznow*vznow);
    vwind   = sqrt((vxnow-xwind)^2 + (vynow-ywind)^2 ...
                + vznow^2);                                  %velocity minus wind component
    reynolds(i) = rhoa(i) * vwind * 2 * rad / visc;               %Reynolds number
    mach(i) = vwind / c_sound;                                    %Mach number
    cdnow(i)=cdnow(i-1);                                          %Cd (using the last value is good enough)

    %Trim arrays
    ifinal = i;
    x=x(1:ifinal); 
    y=y(1:ifinal); 
    z=z(1:ifinal); 
    vx=vx(1:ifinal);
    vy=vy(1:ifinal);
    vz=vz(1:ifinal);
    time=time(1:ifinal);
    mach=mach(1:ifinal);
    cdnow=cdnow(1:ifinal);
    temp=temp(1:ifinal);
    rhoa=rhoa(1:ifinal);
    reynolds=reynolds(1:ifinal);

    %% WRITE OUT RESULT

    if strcmp(write_outfile,'yes')
        %Get current date and time for filename
        now = datestr(datetime);
        filename = sprintf('%s_run%03d.txt',filename_stem,inow);

        fid = fopen(filename, 'w');

        fprintf(fid,'Output from Eject, Matlab version, run on %s local time\n\n',now);
        fprintf(fid,'Block parameters:\n');
        fprintf(fid,'                       Block diameter, m:  %5.2f\n',diam);
        fprintf(fid,'                    Block density, kg/m3:   %4.0f\n',rhor);
        fprintf(fid,'                             Block shape:  %s\n',dragtype);
        fprintf(fid,'                   Initial velocity, m/s: %6.1f\n',vi);
        fprintf(fid,' Ejection angle from horizontal, degrees:   %4.1f\n',thetadeg);
        fprintf(fid,'      Ejection direction, degrees E of N:  %5.1f\n',phideg);
        fprintf(fid,'                         Wind speed, m/s:   %4.1f\n',wind);
        fprintf(fid,'             Wind direction, deg. E of N:   %4.1f\n',winddir_deg);
        fprintf(fid,'landing point meters below takeoff point:   %4.0f\n',xi);
        fprintf(fid,'     extent of reduced drag zone, meters:   %4.1f\n\n',dragred);

        fprintf(fid,'Meteorologic conditions:\n');
        fprintf(fid,'            Temperature at sea level, Celsius:  %4.1f\n',TzeroC);
        fprintf(fid,'            Thermal lapse rate, deg. C per km: %5.2f\n',lapse);
        fprintf(fid,'Elevation of takeoff point above sea level, m:  %4.0f\n',elev);
        fprintf(fid,'                               vent longitude:  %10.5f\n',vent_lon);
        fprintf(fid,'                                vent latitude:  %8.5f\n\n',vent_lat);

        fprintf(fid,'TRAJECTORY CALCULATIONS\n');
        fprintf(fid,'##################################################################################################################################\n');
        fprintf(fid,'    i   t (s)     x (m)     y (m)     z (m)   vx(m/s)  vy (m/s)  vz (m/s)   Mach #     Reynolds #     Cd   Temp (K) rho_air (kg/m3)\n');

        for i=1:ifinal
            fprintf(fid,'%5d  %6.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %8.3f %14.3f %6.3f %10.2f %10.3f\n', ...
               i,time(i),x(i),y(i),z(i),vx(i),vy(i),vz(i),mach(i),reynolds(i),cdnow(i),temp(i),rhoa(i));
        end

        fprintf(fid,'##################################################################################################################################\n');
        fprintf(fid,'Range=%7.1f meters, flight time=%7.1f seconds, final location (x,y,z)=(%7.1f,%7.1f,%7.1f)\n',distance,tfinal,xfinal,yfinal,zfinal);
        fclose(fid);
    end

end

