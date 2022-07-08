%script for computing float energy as a function of float radius

%parallelization
pc = parcluster('local')
parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')))

%number of displacement points
d_num = 8;
%number of radius points
r_num = 9;

options = optimset('PlotFcns',@optimplotx,'MaxIter',100);
float_min = zeros(r_num+1,d_num+1, 2);
energy_min = zeros(r_num+1,d_num+1, 1);

parfor r = 0:r_num

    %hard code in the actual values
    R = .03*(r+1)
    c_stick = 1/tan(50/180*pi);
    ell = 1;
    F = -(2*pi*c_stick+1)*R*R;
    
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
            d = i*(.4/d_num);
            E_min_init = [A*ell^2,-psi_coeff*d];

            f = @(x)radius_sweep_func(x,R,F,c_stick,d,ell);
            [E_min,E_val] = fminsearch(f,E_min_init,options);
            temp_e(i+1) = E_val;
            temp_min(i+1,:) = E_min;

        catch ME
        end
    end

    energy_min(r+1,:) = temp_e;
    float_min(r+1,:,:) = temp_min;

end

save('radius_sweep_min.mat','float_min')
save('radius_sweep_e.mat','energy_min')
