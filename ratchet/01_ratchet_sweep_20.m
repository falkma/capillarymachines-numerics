%script for generating energy landscape of a float in a ratchet channel
%with contact angle 20deg

%parallelization
pc = parcluster('local')
parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')))

%set how many points in the grid for the angle describing the float
turn_num = 269;
%set how many points in the grid for the angle describing the channel
beta_num = 119;
%sets max number of optimization rounds
options = optimset('MaxIter',360);
float_min = zeros(turn_num+1,beta_num+1,2);
%set the contact angle at the channel
c_stick =  1/tan(20/180*pi);

parfor i = 0:turn_num
    
    turn = i/(turn_num+1)*pi;
    E_min_init = [1.];
    temp = zeros(beta_num+1,2);
    
    for j = 0:beta_num

        %minimization of energy over float height for fixed channel/float geometry
        beta = j/(beta_num+1)*pi;
        f = @(t)rectangular_ratchet(t,turn,beta,c_stick);
        [E_min,E_val] = fminsearch(f,E_min_init,options);
        temp(j+1,:) = [E_min,E_val];
        E_min_init(1) = E_min(1);

    end

    float_min(i+1,:,:) = temp;

end

save('ratchet_sweep_20.mat','float_min')
