%function for computing float overlaps

function o_laps = overlap_list(E,radii)

    %assign positions of three floats to individual variables
    x1 = E(1);
    y1 = E(2);
    x2 = E(4);
    y2 = E(5);
    x3 = E(7);
    y3 = E(8);

    %compute pairwise differences
    float_positions = [x1 y1; x2 y2; x3 y3];
    seps = pdist2(float_positions,float_positions);
    o_laps = zeros(1,3);
    
    %compute overlaps assuming that there is only one touching pairwise
    %allow for a gap of .015
    if seps(1,2) < radii(1)+radii(2)+.015
        o_laps(1) = 1;
    else
    end
    
    if seps(1,3) < radii(1)+radii(3)+.015
        o_laps(3) = 1;
    else
    end
    
    if seps(2,3) < radii(2)+radii(3)+.015
        o_laps(2) = 1;
    else
    end
    

end
