function rknow = rk3d(x,y,z,vx,vy,vz,xwind,ywind,dt,rhoa,cdnow,area,mass,grav,rhor)

%Function that uses a Runga-Kutte ODE integrator to integrate the ballistic
%trajectory
%   x = current x location, m
%   y = current y location, m
%   z = current z location, m
%   vx = current x velocity component, m/s
%   vy = current y velocity component, m/s
%   z  = current z velocity component, m/s
%   xwind = x component of wind speed, m/s
%   ywind = y component of wind speed, m/s
%   dt    = time step, s
%   rhoa  = air density, kg/m3
%   cdnow = drag coefficient
%   area  = block frontal area, m2
%   mass  = block mass, kg
%   grav  = gravitational acceleration, m/s2
%   rhor  = block density, kg/m3

    vwind  = sqrt((vx - xwind).^2 + (vy-ywind).^2 + vz.*vz); %horiz. velocity, minus wind
    vvxwd  = vx - xwind;          %x velocity, minus tailwind
    vxwind = vx - xwind;          %x velocity, minus tailwind
    vvywd  = vy - ywind;
    vywind = vy - ywind;
    vvz = vz;

    derivs1 = derivs(vwind,vvxwd,vvywd,vvz,xwind,ywind,dt,rhoa,cdnow,area,mass,grav,rhor);

    vwind = sqrt((vx+0.5.*derivs1(4)).^2 + (vy+0.5.*derivs1(5)).^2 + (vz+0.5.*derivs1(6)).^2);
    vvxwd = vx + 0.5 * derivs1(4);
    vvywd = vy + 0.5 * derivs1(5);
    vvz   = vz + 0.5 * derivs1(6);

    derivs2 = derivs(vwind,vvxwd,vvywd,vvz,xwind,ywind,dt,rhoa,cdnow,area,mass,grav,rhor);       

    vwind = sqrt((vx+0.5*derivs2(4)).^2 + (vy+0.5.*derivs2(6)).^2 + (vz+0.5.*derivs2(6)).^2);
    vvxwd = vx + 0.5 * derivs2(4);
    vvywd = vy + 0.5 * derivs2(5);
    vvz   = vz + 0.5 * derivs2(6);

    derivs3 = derivs(vwind,vvxwd,vvywd,vvz,xwind,ywind,dt,rhoa,cdnow,area,mass,grav,rhor);       

    vwind = sqrt((vx+0.5.*derivs3(4)).^2 + (vy+0.5.*derivs3(5)).^2 + (vz+0.5.*derivs3(6)).^2);
    vvxwd = vx + 0.5 * derivs3(4);
    vvywd = vy + 0.5 * derivs3(5);
    vvz   = vz + 0.5 * derivs3(6);

    derivs4 = derivs(vwind,vvxwd,vvywd,vvz,xwind,ywind,dt,rhoa,cdnow,area,mass,grav,rhor);       

    r = zeros(6,4);
    for i=1:6
        r(i,1) = derivs1(i);
        r(i,2) = derivs2(i);
        r(i,3) = derivs3(i);
        r(i,4) = derivs4(i);
    end

    d=zeros(6,1);
    rknow = zeros(6,1);

    for i=1:6
        d(i) = r(i, 1) / 6 + r(i, 2) / 3 + r(i, 3) / 3 + r(i, 4) / 6;
        rknow(1) = x + d(1);                %new x
        rknow(2) = y + d(2);                %new x
        rknow(3) = z + d(3);                %new z
        rknow(4) = vxwind + d(4) + xwind;   %new vx
        rknow(5) = vywind + d(5) + ywind;   %new vy
        rknow(6) = vz + d(6);               %new vz
    end

    return