%function for computing energy of a set of three floats in a taichi device cross-section

function E = taichi_radial(x,t,container,masses,radii,c,L,mode)

    warning('off','all')
    %set the fixed constants
    %inverse capillary length squared
    invl2 = 1/L/L;
    %mass of float divided by fluid density
    m1 = masses(1);
    m2 = masses(2);
    m3 = masses(3);
    %boundary condition for container
    stick = c;
    %radius of float
    a1=radii(1);
    a2=radii(2);
    a3=radii(3);

    %assign the variable parts of the geometry
    %rotation angle around z-axis centered on float
    theta1 = t(1);
    %rotation angle around x-axis centered on float
    alpha1 = t(2);
    %rotation angle around y-axis centered on float
    gamma1 = t(3);
    %center of float
    x1 = x(1);
    y1 = x(2);
    z1 = x(3);
    
    %rotation angle around z-axis centered on float
    theta2 = t(4);
    %rotation angle around x-axis centered on float
    alpha2 = t(5);
    %rotation angle around y-axis centered on float
    gamma2 = t(6);
    %center of float
    if mode > 0
        sep = a1+a2+.003;
        x2 = x1+sep*cos(x(4));
        y2 = y1+sep*sin(x(4));
    else
        x2 = x(4);
        y2 = x(5);
    end
    z2 = x(6);

    %rotation angle around z-axis centered on float
    theta3 = t(7);
    %rotation angle around x-axis centered on float
    alpha3 = t(8);
    %rotation angle around y-axis centered on float
    gamma3 = t(9);
    %center of float

    %switching from cartesian to radial coordinates if there is a contact
    if mode > 1
        sep = a2+a3+.003;
        x3 = x2+sep*cos(x(7));
        y3 = y2+sep*sin(x(7));
    else
        x3 = x(7);
        y3 = x(8);
    end
    z3 = x(9);

    R1 = rot_mat(alpha1,theta1,gamma1);
    R2 = rot_mat(alpha2,theta2,gamma2);
    R3 = rot_mat(alpha3,theta3,gamma3);
    
    %create and rotate the float
    float_res = 50;
    float_pts = (1:float_res)*2*pi/float_res;
    flat_float1 = [a1*cos(float_pts); a1*sin(float_pts); zeros(1,float_res)];
    flat_float2 = [a2*cos(float_pts); a2*sin(float_pts); zeros(1,float_res)];
    flat_float3 = [a3*cos(float_pts); a3*sin(float_pts); zeros(1,float_res)];
    
    rot_float1 = R1*flat_float1;
    rot_float2 = R2*flat_float2;
    rot_float3 = R3*flat_float3;
    
    p1 = polyshape(rot_float1(1,:)+x1,rot_float1(2,:)+y1);
    p2 = polyshape(rot_float2(1,:)+x2,rot_float2(2,:)+y2);
    p3 = polyshape(rot_float3(1,:)+x3,rot_float3(2,:)+y3);
    
    enc1 = numel(container(inpolygon(rot_float1(1,:)+x1,rot_float1(2,:)+y1,container(1,:),container(2,:))));
    enc2 = numel(container(inpolygon(rot_float2(1,:)+x2,rot_float2(2,:)+y2,container(1,:),container(2,:))));
    enc3 = numel(container(inpolygon(rot_float3(1,:)+x3,rot_float3(2,:)+y3,container(1,:),container(2,:))));
      
    overlap_mat = overlaps([p1 p2 p3]);
    float_overlap = sum(sum(overlap_mat))-trace(overlap_mat);
    
    if float_overlap > 0
        
        E = NaN;
        return
        
    else

        %set up the dirichlet boundary condition for the float
        n1=R1*[0 0 1]';
        fun1 = @(location,state)z1-1/n1(3)*((location.x-x1)*n1(1)+(location.y-y1)*n1(2));
        n2=R2*[0 0 1]';
        fun2 = @(location,state)z2-1/n2(3)*((location.x-x2)*n2(1)+(location.y-y2)*n2(2));
        n3=R3*[0 0 1]';
        fun3 = @(location,state)z3-1/n3(3)*((location.x-x3)*n3(1)+(location.y-y3)*n3(2));

        %projected area for the float
        Area1 = polyarea(rot_float1(1,:),rot_float1(2,:));
        Area2 = polyarea(rot_float2(1,:),rot_float2(2,:));
        Area3 = polyarea(rot_float3(1,:),rot_float3(2,:));

        %actually creating geometries
        full_xlist = [container(1,:) NaN rot_float1(1,:)+x1 NaN rot_float2(1,:)+x2 NaN rot_float3(1,:)+x3];
        full_ylist = [container(2,:) NaN rot_float1(2,:)+y1 NaN rot_float2(2,:)+y2 NaN rot_float3(2,:)+y3];
        pg = polyshape(full_xlist,full_ylist);

        tr = triangulation(pg);
        model = createpde();
        model.geometryFromMesh(tr.Points', tr.ConnectivityList');
        try
            mesh = generateMesh(model,'Hmax',0.07);

            %need to distinguish outer and inner edges somehow, which here we do
            %via distance cutoff
            r_border1=a1*1.01;
            r_border2=a2*1.01;
            r_border3=a3*1.01;

            %setting up boundary conditions by distinguishing inner and outer edges
            %as well as distinguishing edges belonging to individual floats
            outer_edges = [];
            inner_edges1 = [];
            inner_edges2 = [];
            inner_edges3 = [];
            for i = 1:model.Geometry.NumEdges
                nodes = findNodes(mesh,'region','Edge',i);

                nodes_x = mesh.Nodes(1,nodes);

                nodes_y = mesh.Nodes(2,nodes);
                x_x1 = nodes_x - x1;
                y_y1 = nodes_y - y1;
                x_x2 = nodes_x - x2;
                y_y2 = nodes_y - y2;
                x_x3 = nodes_x - x3;
                y_y3 = nodes_y - y3;
                dist1 = mean(sqrt(x_x1.*x_x1+y_y1.*y_y1));

                dist2 = mean(sqrt(x_x2.*x_x2+y_y2.*y_y2));
                dist3 = mean(sqrt(x_x3.*x_x3+y_y3.*y_y3));

                if dist1 < r_border1
                    inner_edges1 = [inner_edges1 i];

                elseif dist2 < r_border2
                    inner_edges2 = [inner_edges2 i];

                elseif dist3 < r_border3
                    inner_edges3 = [inner_edges3 i];

                else
                    outer_edges = [outer_edges i];
                end
            end

            %Apply Boundary Conditions
            applyBoundaryCondition(model,'neumann','Edge',outer_edges,'q',0,'g',stick);
            applyBoundaryCondition(model,'dirichlet','Edge',inner_edges1,'h',1,'r',fun1);
            applyBoundaryCondition(model,'dirichlet','Edge',inner_edges2,'h',1,'r',fun2);
            applyBoundaryCondition(model,'dirichlet','Edge',inner_edges3,'h',1,'r',fun3);

            %Coefficient of the PDE for h
            specifyCoefficients(model,'m',0,'d',0,'c',1,'a',invl2,'f',0);
            h=model.solvepde;

            %acquire FEM matrcies so we can easily compute energy
            FEM = assembleFEMatrices(model);
            U = h.NodalSolution;
            E = .5*(U.'*FEM.K*U + U.'*FEM.A*U) - U.'*FEM.G + (m1*z1+m2*z2+m3*z3+.5*(Area1*z1*z1+Area2*z2*z2+Area3*z3*z3))*invl2;
            
        catch
            E = NaN;
        end
           
    end
    
    function R = rot_mat(a,b,c)
        %compute rotation matrices
        %rotation around x-axis
        Rx=[1 0 0;0 cos(a) -sin(a) ; 0 sin(a) cos(a)];
        %rotation around y-axis
        Ry=[cos(c) 0 sin(c);0 1 0 ; -sin(c) 0 cos(c)];
        %rotation around z-axis
        Rz=[cos(b) -sin(b) 0; sin(b) cos(b) 0; 0 0 1];
        R=Rx*Ry*Rz;
    end

end
