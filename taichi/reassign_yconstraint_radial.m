%function for computing cartesian description to polar description

function init = reassign_yconstraint_radial(E)

    %loads in contacting floats, computes separation between them
    x1 = E(1);
    y1 = E(2);
    x2 = E(4);
    y2 = E(5);
    sep = sqrt((x1-x2)^2+(y1-y2)^2);

    %reassigns variables to be in polar coordinates
    E(4) = asin((y2-y1)/sep);
    E(5) = -1;
    
    init = E;
        
end
