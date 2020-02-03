% Script that makes inputs that will be read by the script
% call_eject3dfun.m

% The inputs include a series of block initial velocities, trajectory
% angles, diameters, and drag shapes.  Each set of inputs is chosen
% randomly from within the following pdfs:
%    vertical angle (thetadeg):  Gaussian pdf with mean=thetadeg_mean
%                                                  sigma=thetadeg_std
%          initial velocity vi:  Gaussian pdf with mean=vi_mean
%                                                 sigma=vi_std
%            for non-vertical trajectories, the initial velocity is 
%            reduced following the relation:
%                            vi = vi * (sin(theta))^2
%    horizontal direction (deg. E of N):  uniform with min=phideg_min
%                                                      max=phideg_max
%    Block diameter:  lognormal with mean=logdiam_mean
%                                   sigma=logdiam=std


%% BLOCK 1:  VALUES THAT DEFINE THE INPUT DISTRIBUTION

%Location of vent
vent_lon = -155.2814;                     %Lat, lon
vent_lat = 19.4080;

%NUMBER OF BLOCKS TO EJECT
n_ejects  = 5;

%OUTPUT FILENAME
filename = 'input\\eject3d_inputs.txt';

%PARAMETERS THAT DEFINE THE RANGE AND DISTRIBUTION OF INITIAL BLOCK SIZES, SHAPES,
%            EJECTION VELOCITIES, AND TRAJECTORY ANGLES

%The values below define the mean and standard
%deviation of a Gaussian distribution of values.

%INITIAL ANGLE FROM HORIZONTAL (Gaussian distribution)
%The resulting values are filtered to ensure that 0 < theta < 90 degrees
thetadeg_mean = 90.;      %Mean ejection angle (measured vertically from horizontal, in degrees)
thetadeg_std  =  40.;      %standard deviation of ejection angle, degrees

%INITIAL VELOCITY (Gaussian distribution)
vi_mean = 250.;                %mean initial velocity (m/s) for vertical blocks
vi_std  = 50.;                %standard deviation of velocity (m/s)

%BLOCK DIAMETER (LOGNORMAL)
logdiam_mean = 0;               %log of mean diameter, meters
logdiam_std  = 0.25;               %log of std. deviation in diameter

%BLOCK DIRECTION (degrees east of north; uniform distribution)
phimin_degrees = 0.;
phimax_degrees = 360.;

%BLOCK TYPE
%The model randomly chooses between the four types listed in the 'dragtype' variable.
dragtypes = {'hicube', 'locube', 'sphere', 'shell'};
%To change the drag types considered, modify the elements of dragtypes.
%For example, to consider only "shell" types, uncomment the line below, and comment out
% the line above.
%dragtypes = {'shell', 'shell', 'shell', 'shell'};

%% BLOCK 3.  CALCULATE RANDOM DISTRIBUTIONS

% PHYSICAL CONSTANTS
pi     =    3.14159;

%Vertical angle
thetadeg  = gaussian(thetadeg_mean,thetadeg_std,n_ejects);          %Normal distribution
for i=1:n_ejects
    thetadeg(i) = max(5.,thetadeg(i));          %set minimum angle to 5 degrees
    if thetadeg(i) > 90.
        thetadeg(i) = 180. - thetadeg(i);        %Make sure everything is between zero and 90
    end
end
theta = pi*thetadeg/180.;                       %convert to radians

%Initial velocity
vi  = gaussian(vi_mean,vi_std,n_ejects);          %Normal distribution
for i=1:n_ejects
    vi(i) =  vi(i) * sin(theta(i))^2;   %account for higher velocities at steeper angles
    vi(i) = max(2.,vi(i));              %make sure nothing has zero or negative ejection velocity
end

%Block diameter (set a lognormal distribution)
logdiam  = gaussian(logdiam_mean,logdiam_std,n_ejects);          %Normal distribution
diam      = 10.^logdiam;

%Initial horizontal direction
phideg  = phimin_degrees + (phimax_degrees-phimin_degrees)*rand(n_ejects,1);

%Block type
typenow   = 4.*rand(n_ejects,1);                 %generate random numbers between 0 and 4
typenow   = 1+floor(typenow);                    %convert them to integers

%% BLOCK 4:  WRITE OUT INPUT FILE

now = datestr(datetime);

fid=fopen(filename,'w');
fprintf(fid,'Inputs for eject3d, written %s\n\n',now);
fprintf(fid,'         vent location:  longitude=%10.5f,  latitude=%7.5f\n',vent_lon,vent_lat);
fprintf(fid,'      number of blocks:  %5d\n\n',n_ejects);
fprintf(fid,'                  VARIABLE            DISTRIBUTION MEAN OR MIN  STD OR MAX\n');
fprintf(fid,' vertical angle, deg. from horizontal   Gaussian      %4.1f       %4.1f\n',thetadeg_mean,thetadeg_std);
fprintf(fid,'     initial velocity, m/s              Gaussian    %6.1f      %5.1f\n',vi_mean,vi_std);
fprintf(fid,'  log10 block diameter, m               Gaussian      %4.1f       %4.1f\n',logdiam_mean,logdiam_std);
fprintf(fid,' initial direction, deg. E of N         Uniform      %5.1f      %5.1f\n\n',phimin_degrees,phimax_degrees);
fprintf(fid,'     block types:  %8s  %8s  %8s  %8s\n\n',dragtypes{1},dragtypes{2},dragtypes{3},dragtypes{4});

fprintf(fid,'INPUTS FOR EACH RUN\n');
fprintf(fid,'***********************************************\n');
fprintf(fid,'   i      vi    diam   theta     phi   dragtype\n');

for i=1:n_ejects
    dragtypenow = dragtypes{typenow(i)};                                   %set dragtypenow
    fprintf(fid,'%4d%8.1f%8.2f%8.1f%8.1f  %8s\n', ...
            i,vi(i),diam(i),thetadeg(i),phideg(i),dragtypenow);
end
fclose(fid);

%% Plot Data

% Plot velocity vectors from vertical
subplot(1,3,1);
for i=1:n_ejects
    xval = vi(i)*cos(pi*thetadeg(i)/180.);
    yval = vi(i)*sin(pi*thetadeg(i)/180.);
    plot([0 xval],[0 yval],'r-'), hold on;
end

thetanow = 0:pi/90:pi/2;
x_envelope = zeros(length(thetanow),1);
y_envelope = zeros(length(thetanow),1);

for j=1:length(thetanow)
    vinow = vi_mean*sin(thetanow(j))^2;
    x_envelope(j) = vinow * cos(thetanow(j));
    y_envelope(j)  = vinow * sin(thetanow(j));
end
plot(x_envelope,y_envelope,'b--');
axis equal;
axis([0 inf 0 inf]);
grid on;
title('Velocity vectors');
xlabel('m/s horizontal');
ylabel('m/s vertical');

% Plot horizontal orientation
subplot(1,3,2)
for i=1:n_ejects
    xval = vi(i)*cos(pi*(90-phideg(i))/180.);
    yval = vi(i)*sin(pi*(90-phideg(i))/180.);
    plot([0 xval],[0 yval],'r-'), hold on;
end
title('horizontal direction');
xlabel('m/s east');
ylabel('m/s north');
grid on
axis equal

% Plot histogram of block diameter
subplot(1,3,3);
histogram(log10(diam));
title('log block diameter');
xlabel('log10(block diameter, m)');
ylabel('number');
grid on
saveas(gcf,'input\eject3d_inputs.jpg','jpg');

fprintf('All done\n');
