This repository contains matlab scripts that generate a series of inputs for 
ballistic block ejection from a volcanic vent, calculate their trajectories, 
and overlay the trajectories over a terrain model of a volcano to determine 
when and where they land.  The code was translated to Matlab from the Visual Basic program Eject! (Mastin, 2001), converted to 3d, and integrated with geographic
information.

The main files are as follows:
MakeInputs.m           script that generates inputs and writes them to input
                           file :
                               input/eject_inputs.jpg
call_eject3dfun.m      script that reads input file, calls a function that
                           calculates trajectories, and writes the locations
                           of the landing points to the file:
                               output\eject3d_output.txt
eject3dfun.m           function that calculates ballistic trajectories and
                           writes the trajectory for each simulation to the file
                                output\Halemauamua_run###.txt,
                           where "###" is the number of each run.
write_kml.m            script that reads output\eject3d_output.txt and writes
                           a kml file showing the landing points of the blocks:
                                output\eject3d_output.txt
write_trajectories.m   script that reads from the output trajectory files and 
                           writes a kml file showing the block trajectories:
                                output\block_trajectories.kml

Utilities and functions:
  deg2utm.m        function that converts degrees lat/lon to utm coordinates. 
                       Written and copyrighted by Rafael Palacios, 2006, but
                       licensed for public use and posted on Matlab File
                       Exchange 
                       https://www.mathworks.com/matlabcentral/fileexchange/10915-deg2utm
                       The copyright notice is include in a comment block
                       within the file.
  derivs.m         called by rk3d.m; calculates first derivatives of position 
                       and velcity along the trajectory.
  drag_hicube.m    function that calculates the drag coefficient of a hi cube 
                       given the Mach number.  A "hi cube" is a cube whose face
                       is oriented in the direction of travel.
  drag_locube.m    function that calculates the drag coefficient of a lo cube 
                       given the Mach number.  A "hi cube" is a cube whose corner
                       is oriented in the direction of travel.
  drag_shell.m    function that calculates the drag coefficient of a G1 artillery
                       shell given the Mach number. 
  drag_sphere.m   function that calculates the drag coefficient of a sphere given 
                       the Reynolds and Mach numbers
  gaussian.m      function that calculates a series of random values following
                       a Gaussian pdf with a specified mean and standard deviation.
  rk3d.m          One-dimensional (ODE) integrator that integrates the equations
                       of motion along the block trajectory.  Called by eject3dfun.m
  tifreader.m     Reads a geotiff file and creates the file terrain.mat, located
                      in the input directory, with geographic data that are 
                      read by the main scripts.  In particular, the geographic 
                      variables written to terrain.mat are as follows:
                           terrain   an NxM matrix of elevations, in meters asl.
                           t_lat     an NxM matrix of latitude values for each point
                                        in terrain
                           t_lon     an NxM matrix of latitude values for each point
                                        in terrain
                           t_x       an NxM matrix giving UTM x values (eastings) for
                                        each point in terrain
                           t_x       an NxM matrix giving UTM y values (northings) for
                                        each point in terrain
                           utmzone   the utm zone.
                       The geotiff used to create data in the example terrain.mat file
                       is located in input\KIL_30m2019\KIL_30M2019.tif.
                       If you wish to use this code on another volcano, you will have to
                       find a geotiff file for that volcano and run this script.  You
                       may also have to modify it depending on the format of the geotiff
                       you are reading.
  utm2deg.m        function that converts utm coordinates to degrees lat/lon.
                       Written and copyrighted by Rafael Palacios, 2006, but
                       licensed for public use and posted on Matlab File
                       Exchange 
                       https://www.mathworks.com/matlabcentral/fileexchange/10915-deg2utm
                       The copyright notice is include in a comment block
                       within the file.
                       
                       
Example input and output files are in the folders input and output                                
                                
To run the model, do the following:
  1)  If you want to keep the old output files, move those files to a new
          directory or rename them.
  2)  Edit variables in block 1 of MakeInputs.m, then run it.
  3)  run call_eject3dfun.m to do the simulations.
  4)  run write_kml.m  to make a kml file of block landing locations.
  5)  run write_trajectories_km.m to make a kml file of block trajectores.
  
  
Reference:
Mastin, L.G., 2001, A simple calculator of ballistic trajectories for blocks 
ejected during volcanic eruptions, U.S. Geological Survey Open-File Report 01-45, 
U.S. Geological Survey, p. 26.

  

                           