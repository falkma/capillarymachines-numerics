%script for refining trajectory segments of floats in the forward stroke of a taichi device

warning('off','all')

pc = parcluster('local')
parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')))

%load in the geometry of the taichi machine channel
load sliced_taichi_interior.mat
%load in data from previously completed forward simulation
load taichi_forward_sim_min.mat

%how many cross-sections of channel to consider
slice_num = 50;
%max number of optimization iterations allolwed
iter_num = 1000;

float_min_copy = float_min;

parfor i = 60+(1:slice_num)
    
    %physical constants
    radii = .2*[1,1,1];
    L = .271;
    c = 1/tan(20/180*pi);
    masses = 0.0354*[1,1,1];

    %sets channel geometry loaded from sliced_taichi_interior.mat
    container = movelist(410-i);
    container = container{1};
    %centers and normalizes units to fit in with physical constant definitions
    container = container - [41.5 9];
    container = [container(:,1)/10,container(:,2)/10].';

    %initialize each frame with the corresponding frame of the already completed forward simulation
    E_min_init = float_min_copy(i,1:9);
    o_laps = overlap_list(E_min_init,radii);

    %goes through code as if there are no overlaps
    %this is important for checking to see if there are
    %constrained floats we need to split up

    %to refine, we allow more time for floats to move without the radial constraint
    if sum(o_laps) > -1
        
        %function to minimize, assuming all tilt angles are 0
        f_zerocon = @(x)taichi_radial(x,[0,0,0,0,0,0,0,0,0],container,masses,radii,c,L,0);
        %start with an initial test run
        options = optimset('PlotFcns',@optimplotx,'MaxIter',100);
        [E_min,E_val] = fminsearch(f_zerocon,E_min_init,options);
        E_min_init = E_min;
        %recompute overlaps
        o_laps = overlap_list(E_min_init,radii);
        disp([string(i) 'initial unconstrained'])
        disp(o_laps)
        disp(E_min_init)

        %if there are still no overlaps, continue as a normal simulation
        if sum(o_laps) < 1
            options = optimset('PlotFcns',@optimplotx,'MaxIter',iter_num-100);
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

    temp_min = E_min_init;
    temp_E = E_val;

    float_min(i,:) = [temp_min,temp_E];

end

save('taichi_forward_sim_min_refine.mat','float_min')
