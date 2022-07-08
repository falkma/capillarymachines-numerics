%function for computing energy of capillary trap with specified geometry and float radius

function E = radius_sweep_func(x,a_t,m,c,center,cap_length)

    %set the fixed constants
    %inverse capillary length squared
    ell = cap_length;
    invl2 = 1/(ell*ell);
    %mass of float divided by fluid density
    F = m;
    %boundary condition for container
    stick = c;
    %radius of container
    r1 = 1;
    %radius of float
    a=a_t;

    %assign the variable parts of the geometry
    %rotation angle around z-axis centered on float
    theta = 0;
    %rotation angle around x-axis centered on float
    alpha = 0;
    %rotation angle around y-axis centered on float
    gamma = x(2);
    %center of float
    x0 = 0;
    y0 = 0;
    z0 = x(1);

    %container geometry
    C_cont = [1 center 0 r1 0 0];

    %compute rotation matrices
    %rotation around x-axis
    Rx=[1 0 0;0 cos(alpha) -sin(alpha) ; 0 sin(alpha) cos(alpha)];
    %rotation around y-axis
    Ry=[cos(gamma) 0 sin(gamma);0 1 0 ; -sin(gamma) 0 cos(gamma)];
    %rotation around z-axis
    Rz=[cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
    R=Rx*Ry*Rz;
    
    %rotated ellipse size projection
    a_rot = R*[a 0 0]';
    b_rot = R*[0 a 0]';
    a1 = sqrt(a_rot(1)^2+a_rot(2)^2);
    b1 = sqrt(b_rot(1)^2+b_rot(2)^2);
    
    %float geometry
    C_float = [4 x0 y0 a1 b1 theta];

    %set up the dirichlet boundary condition for the float
    n=R*[0 0 1]';
    fun= @(location,state)z0-1/n(3)*((location.x-x0)*n(1)+(location.y-y0)*n(2));

    %projected area for the float
    Area=pi*a1*b1;

    try
        % combining geometries
        gd = [C_cont; C_float]';
        % Names for the two geometric objects
        ns = (char('C_cont','C_float'))';
        % Set formula
        sf = 'C_cont-C_float';

        %actually creating geometries
        [dl, bt]=decsg(gd,sf,ns);
        [dl2,~] = csgdel(dl,bt);
        model=createpde();
        geometryFromEdges(model,dl2);

        %need to distinguish outer and inner edges somehow, which here we do
        %via distance cutoff
        r2=a_t*1.05;
        %r2=r1*.95;
        inner_edges = find((abs(dl2(4,:))).^2+(abs(dl2(2,:))).^2 <r2^2);
        outer_edges = find((abs(dl2(4,:))).^2+(abs(dl2(2,:))).^2 >=r2^2);

        %Apply Boundary Conditions 
        applyBoundaryCondition(model,'neumann','Edge',outer_edges,'q',0,'g',stick);
        %applyBoundaryCondition(model,'neumann','Edge',inner_edges,'q',0,'g',-stick);
        applyBoundaryCondition(model,'dirichlet','Edge',inner_edges,'h',1,'r',fun);

        %Coefficient of the PDE for h
        specifyCoefficients(model,'m',0,'d',0,'c',1,'a',invl2,'f',0);
        generateMesh(model,'Hmax',0.07);
        h=model.solvepde;
        FEM = assembleFEMatrices(model);
        U = h.NodalSolution;
        %E = .5*(U.'*FEM.K*U + U.'*FEM.A*U) - U.'*FEM.G + mdivrho*z0*invl2 + .5*Area*z0*z0*invl2 + 3.14/8*invl2*a*a*a*a*sin(gamma)*sin(gamma);
        E = .5*(U.'*FEM.K*U + U.'*FEM.A*U) - U.'*FEM.G - F*z0 + .5*Area*z0*z0*invl2 + pi/8*invl2*a*a*a*a*gamma*gamma;
    
        figure(2)
        pdegplot(model)
        figure(3)
        pdeplot(model,'ZData',U);
        
    catch
        E = NaN;
    end

end
