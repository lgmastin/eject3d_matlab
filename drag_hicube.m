function drag = drag_hicube(mach)

    % Drag of a hi cube (i.e. a cube whose face is facing toward the
    % front).  Drag function as a function of the Mach number, with the
    % relationship digitized from Fig. 17 (p. 16-14) or
    %Hoerner, S.F., 1965, Fluid Dynamic Drag: Vancouver, WA, 
    %  published by the Author.
    % input = Mach number
    % Output = drag coefficient at that Mach number

    hicube_mach = [-0.9624766, -0.8649155, -0.6435272, -0.5159475, ...
                            -0.3921201, -0.2645403, -0.1932458, -0.1181989, ...
                            -0.07317073,-0.03564728, 0.001876173,0.0619137, ...
                            0.1332083,   0.2345216,  0.3245779,  0.4333959, ...
                            0.4934334,   0.5234522,  0.9362102];
    hicube_cd  = [1.064151,    1.056604,   1.064151,   1.086792,  ...
                           1.116981,    1.184906,   1.237736,   1.320755,  ...
                           1.396226,    1.486792,   1.607547,   1.728302,  ...
                           1.788679,    1.811321,   1.811321,   1.773585,  ...
                           1.743396,    1.728302,   1.720755];
    hicube_mach = 10.^hicube_mach;
    %Create interpolation function
    if mach < hicube_mach(1)
        drag = hicube_cd(1);
    elseif mach > hicube_mach(19)
        drag = hicube_cd(19);
    else
        %Create interpolation function
        drag = interp1(hicube_mach,hicube_cd,mach,'pchip');
    end
end
    
    
