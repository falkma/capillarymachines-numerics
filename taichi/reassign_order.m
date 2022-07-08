%function for reassigning float indices so overlapping floats are first

function [init_E,init_r,init_m] = reassign_order(E,radii,masses,o_laps)

    %shift ordering of floats so that any overlaps are first
    perm_num = 0;
    while o_laps(1) ~= 1
        o_laps = circshift(o_laps,1);
        perm_num = perm_num + 1;
    end

    %shift other variables according to this order
    init_E = circshift(E,3*perm_num);
    init_r = circshift(radii,perm_num);
    init_m = circshift(masses,perm_num);
end
