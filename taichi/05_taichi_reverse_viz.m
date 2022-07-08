%script to load and visualize the trajectories of floats in the reverse taichi simulation

load sliced_taichi_interior.mat
load taichi_reverse_sim_min.mat
slice_num = 400;

for i = 0+(1:slice_num)
    container = movelist(10+i);
    container = container{1};
    container = container - [41.5 9];
    container = [container(:,1)/10,container(:,2)/10].';
    radii = .2*[1,1,1];
    L = 0.271;
    c = 1/tan(20/180*pi);
    masses = 0.0354*[1,1,1];
    taichi_radial(float_min(i,1:9),[0,0,0,0,0,0,0,0,0],container,masses,radii,c,L,0);
    filename = sprintf('taichi_reverse_movie/taichi_movie%d.png',i);
    saveas(gcf,filename)
end