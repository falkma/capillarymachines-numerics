%script for computing float energy as a function of channel contact angle

%parallelization
pc = parcluster('local')
parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')))

%number of displacement points
d_num = 8;
%number of slope points
c_num = 16;

options = optimset('PlotFcns',@optimplotx,'MaxIter',100);
float_min = zeros(c_num+1,d_num+1, 2);
energy_min = zeros(c_num+1,d_num+1, 1);

parfor c = 0:c_num

    %hard code in the actual values
    R = .1
    angle = 10 + 5*c
    c_stick = 1/tan(angle/180*pi);
    ell = 1;
    F = -.2;
    
    temp_min = zeros(d_num+1, 2);
    temp_e = zeros(d_num+1, 1);

    %initialize based on theory
    B0 = -F/2/pi;
    A = 2*(c_stick - B0);
    A1 = (2*B0 - c_stick)/(1-R*R);
    B1 = R*R*A1;
    psi_coeff = A1 + B1/R/R;

    for i = 0:d_num

        try
            d = i*(.1/d_num);
            E_min_init = [A*ell^2,-psi_coeff*d];

            f = @(x)slope_sweep_func(x,R,F,c_stick,d,ell);
            [E_min,E_val] = fminsearch(f,E_min_init,options);
            temp_e(i+1) = E_val;
            temp_min(i+1,:) = E_min;

        catch ME
        end
    end

    energy_min(c+1,:) = temp_e;
    float_min(c+1,:,:) = temp_min;

end

save('slope_sweep_min.mat','float_min')
save('slope_sweep_e.mat','energy_min')
