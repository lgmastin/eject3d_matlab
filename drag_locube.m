function drag = drag_locube(mach)

    % Drag of a lo cube (i.e. a cube whose corner is facing toward the
    % front).  Drag function as a function of the Mach number, with the
    % relationship digitized from Fig. 17 (p. 16-14) or
    %Hoerner, S.F., 1965, Fluid Dynamic Drag: Vancouver, WA, 
    %  published by the Author.
    % input = Mach number
    % Output = drag coefficient at that Mach number

    %Digitized curve of drag coefficient versus Mach number
    locube_mach= [-0.9623352, -0.7702448, -0.6007533, -0.4387947, ...
                  -0.3182674, -0.2165725, -0.1412429, -0.08097929, ...
                  -0.03954802, 0.01318267, 0.06214689, 0.126177, ...
                   0.1939736,   0.2730697,  0.3559322, 0.4274953, ...
                   0.5028248, 0.5969868, 0.9397364];
    locube_cd = [0.7698113, 0.7698113, 0.7773585, 0.7849057, ...
                  0.8150944,   0.8528302,  0.8981132,   0.9358491,  ...
                  0.9962264,   1.109434,   1.184906,    1.222641,   ...
                  1.237736,    1.230189,   1.207547,    1.177359,   ...
                  1.14717,     1.139623,   1.139623];
    locube_mach = 10.^locube_mach;
    if mach < locube_mach(1)
        drag = locube_cd(1);
    elseif mach > locube_mach(19)
        drag = locube_cd(19);
    else
        %interpolate using the given mach number
        drag = interp1(locube_mach,locube_cd,mach,'pchip');
    end
end
  
