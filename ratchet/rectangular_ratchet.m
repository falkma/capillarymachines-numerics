%function for computing energy of float configurations in a ratchet geometry

function E = rectangular_ratchet(t,angle,beta,c_stick)

    %set the fixed constants
    %inverse capillary length squared
    invl2 = 13.61;
    %mass of float divided by fluid density
    mdivrho = 0.02;
    %boundary condition for container
    stick = c_stick;
    %dimensions of rectangular container
    r1 = .4;
    r2 = .57;
    %half-widths of rectangular channel through the container
    x_o = 1.0;
    y_o = .1915;
    %half-widths of float
    x_i = .25;
    y_i = .15;
    
    %assign the variable parts of the geometry
    %rotation angle around z-axis centered on float
    theta = angle;
    %rotation angle around x-axis centered on float
    alpha = 0;
    %rotation angle around y-axis centered on float
    gamma = 0;
    %center of float
    x0 = 0;
    y0 = 0;
    z0 = t(1);
    %rotation angle around z-axis centered on container ellipse
    %beta = x(6);

    %rectangular container geometry
    %C_ellipse = [4 0 0 r2 r1 beta 0 0 0 0];
    init_rectangle = [-r2 r2 r2 -r2; -r1 -r1 r1 r1];
    rot1 = [cos(beta) -sin(beta) ; sin(beta) cos(beta)];
    rot_C2 = rot1*init_rectangle;
    C_rect=[3 4 rot_C2(1,:) rot_C2(2,:)];
    
    %rectangular channel geometry
    C_channel = [3 4 -x_o x_o x_o -x_o -y_o -y_o y_o y_o];
    
    %compute rotation matrices
    %rotation around x-axis
    Rx=[1 0 0;0 cos(alpha) -sin(alpha) ; 0 sin(alpha) cos(alpha)];
    %rotation around y-axis
    Ry=[cos(gamma) 0 sin(gamma);0 1 0 ; -sin(gamma) 0 cos(gamma)];
    %rotation around z-axis
    Rz=[cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
    R=Rx*Ry*Rz;
    
    %rotating the float
    flat_float = [-x_i x_i x_i -x_i; -y_i -y_i y_i y_i; 0 0 0 0];
    rot_float = R*flat_float;
    
    %float geometry
    %we project out the z-dimension
    C_float = [3 4 rot_float(1,:)+x0 rot_float(2,:)+y0];

    %set up the dirichlet boundary condition for the float
    n = R*[0 0 1]';
    fun = @(location,state)z0-1/n(3)*((location.x-x0)*n(1)+(location.y-y0)*n(2));

    %projected area for the float
    Area = polyarea(rot_float(1,:),rot_float(2,:));

    % combining geometries
    gd = [C_rect; C_channel; C_float]';
    % Names for the two geometric objects
    ns = (char('C_rect','C_channel','C_float'))';
    % Set formula
    sf = '(C_rect+C_channel)-C_float';

    %actually creating geometries
    [dl, bt]=decsg(gd,sf,ns);
    [dl2,~] = csgdel(dl,bt);
    model=createpde();
    geometryFromEdges(model,dl2);

    %need to distinguish outer and inner edges somehow, which here we do
    %via distance cutoff
    r_border=r1*.98;
    inner_edges = find((abs(dl2(4,:))).^2+(abs(dl2(2,:))).^2 <r_border^2);
    outer_edges = find((abs(dl2(4,:))).^2+(abs(dl2(2,:))).^2 >=r_border^2);

    %Apply Boundary Conditions 
    applyBoundaryCondition(model,'neumann','Edge',outer_edges,'q',0,'g',stick);
    applyBoundaryCondition(model,'dirichlet','Edge',inner_edges,'h',1,'r',fun);

    %Coefficient of the PDE for h
    specifyCoefficients(model,'m',0,'d',0,'c',1,'a',invl2,'f',0);
    generateMesh(model,'Hmax',0.07)
    h=model.solvepde;
    FEM = assembleFEMatrices(model);
    U = h.NodalSolution;
    E = .5*(U.'*FEM.K*U + U.'*FEM.A*U) - U.'*FEM.G + mdivrho*z0*invl2 + .5*Area*z0*z0*invl2;
end
