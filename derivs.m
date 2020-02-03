    function derivsnow = derivs( vwind,vvxwd,vvywd,vvz,xwind,ywind,dt,rhoa,cdnow,area,mass,grav,rhor )
    
    %  Function called by subroutine rk3d, which calculates the first
    %  derivaties of x, y, z, dx, dy, and dz.  Inputs are:
    %    vwind =  horiz. velocity, minus wind, m/s
    %    vvxwd =  x component of horizontal velocity, minus wind, m/s
    %    vvywd =  y component . .. 
    %    vvz   =  x component of velocity (wind is assumed to be only horizontal)
    %    xwind =  x component of wind, m/s
    %    ywind =  y component of wind, m/s
    %       dt =  time step (s)
    %     rhoa =  air density, kg/m3
    %    cdnow =  block drag coefficient
    %     area =  block frontal area, m2
    %     mass =  block mass, kg
    %     grav =  gravitation accelecation (9.81 m/s2)
    %     rhor =  block density, kg/m3
 
    derivsnow = zeros(6,1);

    derivsnow(1) = (vvxwd + xwind) * dt;                                                 %distance travelled in x
    derivsnow(2) = (vvywd + ywind) * dt;                                                 %distance travelled in y
    derivsnow(3) = vvz * dt;                                                             %distance travelled in z
    derivsnow(4) = -rhoa * cdnow * area * vwind ^ 2 * vvxwd / (2 * mass * vwind) * dt;  %change in vx
    derivsnow(5) = -rhoa * cdnow * area * vwind ^ 2 * vvywd / (2 * mass * vwind) * dt;  %change in vy
    derivsnow(6) = (-grav * (rhor - rhoa) / rhor - rhoa * cdnow * area * vwind ^ 2 * vvz /(2 * mass * vwind)) * dt;  %in vz

    return
