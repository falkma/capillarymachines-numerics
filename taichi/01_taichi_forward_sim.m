%script for computing trajectory of floats in the forward stroke of a taichi device

warning('off','all')

%load in the geometry of the taichi machine channel
load sliced_taichi_interior.mat

%how many cross-sections of channel to consider
slice_num = 400;
%max number of optimization iterations allolwed
iter_num = 200;

%physical constants
radii = .2*[1,1,1];
L = .271;
c = 1/tan(20/180*pi);
masses = 0.0354*[1,1,1];

%initialize simulation with a guess of the last positions of the reverse simulation
E_min_init = [-1.3728 0.0085 -0.0375 -0.9698 0.0085 -0.0374 2.4588 0.0085 -0.0190];
o_laps = overlap_list(E_min_init,radii);
float_min = zeros(slice_num,10);

for i = 0+(1:slice_num)
    
    %sets channel geometry loaded from sliced_taichi_interior.mat
    container = movelist(410-i);
    container = container{1};
    %centers and normalizes units to fit in with physical constant definitions
    container = container - [41.5 9];
    container = [container(:,1)/10,container(:,2)/10].';

    %goes through code as if there are no overlaps
    %this is important for checking to see if there are
    %constrained floats we need to split up
    if sum(o_laps) > -1
        
        %function to minimize, assuming all tilt angles are 0
        f_zerocon = @(x)taichi_radial(x,[0,0,0,0,0,0,0,0,0],container,masses,radii,c,L,0);
        %start with an initial test run
        options = optimset('PlotFcns',@optimplotx,'MaxIter',10);
        [E_min,E_val] = fminsearch(f_zerocon,E_min_init,options);
        E_min_init = E_min;
        %recompute overlaps
        o_laps = overlap_list(E_min_init,radii);
        disp([string(i) 'initial unconstrained'])
        disp(o_laps)
        disp(E_min_init)

        %if there are still no overlaps, continue as a normal simulation
        if sum(o_laps) < 1
            options = optimset('PlotFcns',@optimplotx,'MaxIter',iter_num-10);
            [E_min,E_val] = fminsearch(f_zerocon,E_min_init,options);
            E_min_init = E_min;
            o_laps = overlap_list(E_min_init,radii);
            disp([string(i) 'final unconstrained'])
            disp(o_laps)
            disp(E_min_init)
        else
        end

    else        
    end
    
    %if there are overlaps, we recompute geometry, and run simulation in constraint mode
    if sum(o_laps) > 0
        
        options = optimset('PlotFcns',@optimplotx,'MaxIter',iter_num);
        %shuffle orders so touching floats are in a particular order
        [E_min_init,radii,masses] = reassign_order(E_min_init,radii,masses,o_laps);
        %now convert the cartesian geometry to a polar geometry
        E_min_constrained = reassign_yconstraint_radial(E_min_init);
        %run the optimization function in constraint mode
        f_onecon = @(x)taichi_radial(x,[0,0,0,0,0,0,0,0,0],container,masses,radii,c,L,1);
        [E_min,E_val] = fminsearch(f_onecon,E_min_constrained,options);
        %convert the polar geometry back to cartesian geometry
        E_min_init = reassign_y_radial(E_min,radii);
        %recompute overlaps
        o_laps = overlap_list(E_min_init,radii);
        disp([string(i) 'constrained'])
        disp(o_laps)
        disp([E_min_init])

    else
    end
        
    %saves geometry to file
    float_min(i,:) = [E_min_init,E_val];
    save('taichi_forward_sim_min.mat','float_min')


end

save('taichi_forward_sim_min.mat','float_min')