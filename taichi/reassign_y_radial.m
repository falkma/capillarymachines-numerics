%function for converting polar description to cartesian description

function init = reassign_y_radial(E,radii)

    %rewrite variables on which optimization is being performed
    x1 = E(1);
    y1 = E(2);
    theta = E(4);

    %allow a small buffer of .003
    sep = radii(1)+radii(2)+.003;

    %reassigns variables to be in cartesian coordinates
    E(4) = x1+sep*cos(theta);
    E(5) = y1+sep*sin(theta);
    init = E;

end
