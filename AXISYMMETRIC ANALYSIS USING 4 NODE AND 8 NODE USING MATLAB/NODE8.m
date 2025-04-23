clf;
clc;
clear;
%% Parameters   
ri = 30; % Inner radius of the cylinder (m)
ro = 50; % Outer radius of the cylinder (m)
L = 40; % Length of the cylinder (m)
E = 2e5;   % Young's modulus (Pa)
v = 0.3;   % Poisson's ratio
P = 40;    % Pressure applied on the right side (N/m²)

% Constitutive matrix
C = (E / ((1+v)*(1-2*v))) * [1-v, v, v, 0; 
                             v, 1-v, v, 0; 
                             v, v, 1-v, 0;
                             0, 0, 0, (1-2*v)/2];

% Define node coordinates
nodes=readmatrix('coords.csv');
nodes(:,3)=nodes(:,3)-ri;

elements1 = readmatrix('element connectivity.csv');

% Extract main nodes (first 4 columns)
elements = elements1(:, 1:4);

% Extract mid-side nodes (last 4 columns)
elements2 = elements1(:, 5:8);

%% Plot the mesh
figure;
hold on;
for i = 1:size(elements,1)
    elemNodes = elements(i, :);
    r = nodes(elemNodes, 2);
    z = nodes(elemNodes, 3);
    fill(r, z, 'c', 'FaceAlpha', 0.3, 'EdgeColor', 'k', 'LineWidth', 1.5);
end
% Plot nodes
scatter(nodes(:,2), nodes(:,3), 50, 'r', 'filled');
for i = 1:size(nodes,1)
    text(nodes(i,2), nodes(i,3), num2str(nodes(i,1)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10);
end

grid on;
axis equal;
xlabel('X Coordinate');
ylabel('Y Coordinate');
title('Mesh Visualization');
hold off;
    %% Define Symbolic Shape Functions for 8-Node Quadrilateral
    syms xi eta
    
    % Shape functions for 8-node serendipity quadrilateral
    N1 = -(1 - xi) * (1 - eta) * (1 + xi + eta) / 4;
    N2 = -(1 + xi) * (1 - eta) * (1 - xi + eta) / 4;
    N3 = -(1 + xi) * (1 + eta) * (1 - xi - eta) / 4;
    N4 = -(1 - xi) * (1 + eta) * (1 + xi - eta) / 4;
    N5 = (1 - xi^2) * (1 - eta) / 2; % Mid-edge bottom
    N6 = (1 + xi) * (1 - eta^2) / 2; % Mid-edge right
    N7 = (1 - xi^2) * (1 + eta) / 2; % Mid-edge top
    N8 = (1 - xi) * (1 - eta^2) / 2; % Mid-edge left
    
    Ni = [N1, N2, N3, N4, N5, N6, N7, N8];
    
    % Compute derivatives symbolically
    dN_dxi = [- ((eta - 1)*(eta + xi + 1))/4 - ((eta - 1)*(xi - 1))/4, 
        ((eta - 1)*(eta - xi + 1))/4 - ((eta - 1)*(xi + 1))/4, 
        ((eta + 1)*(eta + xi - 1))/4 + ((eta + 1)*(xi + 1))/4, 
        ((eta + 1)*(xi - eta + 1))/4 + ((eta + 1)*(xi - 1))/4, 
        xi*(eta - 1), 
        1/2 - eta^2/2, 
        -xi*(eta + 1), 
        eta^2/2 - 1/2];   
    dN_deta = [- ((xi - 1)*(eta + xi + 1))/4 - ((eta - 1)*(xi - 1))/4, 
        ((xi + 1)*(eta - xi + 1))/4 + ((eta - 1)*(xi + 1))/4, 
        ((xi + 1)*(eta + xi - 1))/4 + ((eta + 1)*(xi + 1))/4, 
        ((xi - 1)*(xi - eta + 1))/4 - ((eta + 1)*(xi - 1))/4,
        xi^2/2 - 1/2, 
        -eta*(xi + 1), 
        1/2 - xi^2/2,
        eta*(xi - 1)]; 
    
     %% Gauss Quadrature Points and Weights (3x3 for 8-Node Elements)
   gauss_points = [-sqrt(3/5), 0, sqrt(3/5)]; % Standard points for 3-point quadrature
weights = [5/9, 8/9, 5/9]; % Weights are all 1 for 3×3 integration
    %% Initialize Global Stiffness Matrix
    num_nodes = size(nodes,1);
    num_dof = 2 * num_nodes;  
    K_global = zeros(num_dof, num_dof);
    F_global = zeros(num_dof, 1);
    num_elements = size(elements,1); % number of elements
    node_coords(:,1)=nodes(:,2); % storing r coords in new variable
    node_coords(:,2)=nodes(:,3); % storing z coords in new variable
   
    %% Compute Stiffness Matrix for Axisymmetric 8-Node Elements
    for elem = 1:num_elements
        node_ids = elements1(elem, :);
        coords = node_coords(node_ids, :);
        Ke = zeros(16, 16);
    
        for i = 1:3  % Use 3x3 Gauss Quadrature
            for j = 1:3
                xi_val = gauss_points(i);
                eta_val = gauss_points(j);
                
                % Evaluate shape function derivatives at Gauss points
                dN_dxi_val = double(subs(dN_dxi, [xi, eta], [xi_val, eta_val]));
                dN_deta_val = double(subs(dN_deta, [xi, eta], [xi_val, eta_val]));
    
                % Compute Jacobian Matrix 
               J = [dN_dxi_val'; dN_deta_val'] * coords;
                detJ = det(J);
                Jinv = inv(J);
                
                % Compute B matrix
                B = zeros(4, 16);
                Ni_val = double(subs(Ni, [xi, eta], [xi_val, eta_val]));
                r_gauss = Ni_val * coords(:,1);
                
                for k = 1:8  %CONSTRUCTION OF B MATRIX
                    dN_nat = [dN_dxi_val(k); dN_deta_val(k)];
                    dN_xy = Jinv * dN_nat;
                    
                    B(1, 2*k-1) = dN_xy(1);  
                    B(2, 2*k)   = dN_xy(2);  
                    B(3, 2*k-1) = Ni_val(k) / r_gauss;  
                    B(4, 2*k-1) = dN_xy(2);  
                    B(4, 2*k)   = dN_xy(1);  
                end
                
                % Stiffness matrix at Gauss point
                Ke = Ke + B' * C * B * detJ * weights(i) * weights(j) * 2 * pi * r_gauss;
            end
        end
        
        % Assemble into Global Stiffness Matrix
        dof_map = [2*node_ids-1; 2*node_ids];
        dof_map = dof_map(:)';    
        K_global(dof_map, dof_map) = K_global(dof_map, dof_map) + Ke;
    end
%% Apply boundary conditions
fixed_dof = [];
for i = 1:size(nodes, 1) % 21 numbeer of nodes
    if (nodes(i, 3) == 0 || nodes(i, 3) == max(nodes(:, 3))) % Fix top and bottom plane
        fixed_dof = [fixed_dof; i*2]; % fixing Uz dof
    end
end

for i=1:size(fixed_dof,1)
K_global(fixed_dof(i), :) = 0;  
K_global(:, fixed_dof(i)) = 0;  

K_global((fixed_dof(i)), (fixed_dof(i))) = 1;
F_global(fixed_dof(i)) = 0; % Ensure force vector consistency
end

%% Apply Internal Pressure Load at Inner Radius (r = ri)
force_dof = [];
for i = 1:size(nodes, 1) % number of nodes
    if (nodes(i, 2) == ri) 
        force_dof = [force_dof; i]; % force applied in r direction only
    end
end
for i = 1:size(nodes, 1)
    if nodes(i, 3)~= 0
        zzz = i + 2;
        break;
    end
end
zzzz=nodes(zzz,3); %element size in z direction
full_force = P * 2 * pi * ri *zzzz;
force_dof2=2*force_dof-1; %applying force only to r direction
for i = 1:size(force_dof2,1)
    if i < 3 % applying at the end nodes
        F_global(force_dof2(i)) = F_global(force_dof2(i)) + full_force*(1/6); % Half force for first two DOFs
    elseif any(elements2(:) == force_dof(i)) % check for mid nodes
        F_global(force_dof2(i)) = F_global(force_dof2(i)) + full_force*(4/6);
    else
        F_global(force_dof2(i)) = F_global(force_dof2(i)) + full_force*(2/6); % connecting nodes
    end
end

U = K_global\F_global; % Solve for displacement

Ur = []; % Initialize an empty array to store Ux values

for i = 1:2:(size(F_global,1)-1) %all odd terms are Ux
    Ur = [Ur; U(i)]; % Append values instead of overwriting
end
Uz=[];
for i = 2:2:size(F_global,1)
    Uz = [Uz; U(i)]; % Append values instead of overwriting
end
maxUx=max(Ur)
minUx=min(Ur)
maxUy=max(Uz)
minUy=min(Uz)


%% 

% Solve for displacement
U = K_global\F_global;
% Extract displacement components
Ur = U(1:2:end); % X-displacement (radial)
Uz = U(2:2:end); % Y-displacement (axial)

% Compute displacement magnitude
U_magnitude = sqrt(Ur.^2 + Uz.^2);

% Plot displacement contour
figure;
hold on;
for i = 1:size(elements,1)
    elemNodes = elements(i, :);
    r = nodes(elemNodes, 2);
    z = nodes(elemNodes, 3);
    c = U_magnitude(elemNodes); % Displacement magnitude at element nodes
    patch(r, z, c, 'EdgeColor', 'k', 'FaceColor', 'interp');
end
scatter(nodes(:,2), nodes(:,3), 50, U_magnitude, 'filled'); % Overlay node colors
colorbar;
title('Displacement Contour');
xlabel('X Coordinate');
ylabel('Y Coordinate');
axis equal;
grid on;
hold off;
%% Initialize Storage
stress_r_nodal = zeros(num_nodes, 1);
stress_z_nodal = zeros(num_nodes, 1);
stress_theta_nodal = zeros(num_nodes, 1);
node_count = zeros(num_nodes, 1);

%% Compute Nodal Stresses
for elem = 1:num_elements
    node_ids = elements1(elem, :);
    coords = node_coords(node_ids, :);
    elem_dof = [2*node_ids-1; 2*node_ids]; 
    elem_dof = elem_dof(:)';  
    Ue = U(elem_dof);  

    for local_node = 1:8  % Loop over element nodes
        xi_val = [-1, 1, 1, -1, 0, 1, 0, -1]; % Natural coords of 8 nodes
        eta_val = [-1, -1, 1, 1, -1, 0, 1, 0];

        % Compute Shape Function Derivatives
        dN_dxi_val = double(subs(dN_dxi, [xi, eta], [xi_val(local_node), eta_val(local_node)]));
        dN_deta_val = double(subs(dN_deta, [xi, eta], [xi_val(local_node), eta_val(local_node)]));

        % Compute Jacobian & B Matrix
        J = [dN_dxi_val'; dN_deta_val'] * coords;
        Jinv = inv(J);
        Ni_val = double(subs(Ni, [xi, eta], [xi_val(local_node), eta_val(local_node)]));
        r_node = Ni_val * coords(:,1); % Compute r at node
        
        B = zeros(4, 16);
        for k = 1:8  
            dN_nat = [dN_dxi_val(k); dN_deta_val(k)];
            dN_xy = Jinv * dN_nat;
            
            B(1, 2*k-1) = dN_xy(1);  % ε_rr
            B(2, 2*k)   = dN_xy(2);  % ε_zz
            B(3, 2*k-1) = Ni_val(k) / r_node;  % ε_θθ
            B(4, 2*k-1) = dN_xy(2);  
            B(4, 2*k)   = dN_xy(1);  
        end

        % Compute Stresses
        strain = B * Ue;
        stress = C * strain;
        
        % Assign stresses to corresponding global node
        global_node = node_ids(local_node);
        stress_r_nodal(global_node) = stress_r_nodal(global_node) + stress(1);
        stress_z_nodal(global_node) = stress_z_nodal(global_node) + stress(2);
        stress_theta_nodal(global_node) = stress_theta_nodal(global_node) + stress(3);
        node_count(global_node) = node_count(global_node) + 1;
    end
end

% Average stress values at shared nodes
stress_r_nodal = stress_r_nodal ./ node_count;
stress_z_nodal = stress_z_nodal ./ node_count;
stress_theta_nodal = stress_theta_nodal ./ node_count;

%% Display Maximum and Minimum Stresses
fprintf('Maximum Radial Stress: %.4f MPa\n', max(stress_r_nodal));
fprintf('Minimum Radial Stress: %.4f MPa\n', min(stress_r_nodal));
fprintf('Maximum Axial Stress: %.4f MPa\n', max(stress_z_nodal));
fprintf('Minimum Axial Stress: %.4f MPa\n', min(stress_z_nodal));
fprintf('Maximum Hoop Stress: %.4f MPa\n', max(stress_theta_nodal));
fprintf('Minimum Hoop Stress: %.4f MPa\n', min(stress_theta_nodal));

%% Extract X, Y coordinates
X = node_coords(:, 1);  
Y = node_coords(:, 2);  

%% Plot Stress Contours Using Patch (Ignoring Middle Nodes)
corner_indices = [1, 2, 3, 4]; % Only use corner nodes (no midside nodes)

figure;
hold on;
colormap jet;
for e = 1:num_elements
    nodes = elements1(e, corner_indices); % Use only corner nodes
    patch(X(nodes), Y(nodes), stress_r_nodal(nodes), 'FaceColor', 'interp', 'EdgeColor', 'none');
end
colorbar;
title('Radial Stress Contour (\sigma_r)');
xlabel('X'); ylabel('Y');
axis equal; hold off;

figure;
hold on;
colormap jet;
for e = 1:num_elements
    nodes = elements1(e, corner_indices);
    patch(X(nodes), Y(nodes), stress_z_nodal(nodes), 'FaceColor', 'interp', 'EdgeColor', 'none');
end
colorbar;
title('Axial Stress Contour (\sigma_z)');
xlabel('X'); ylabel('Y');
axis equal; hold off;

figure;
hold on;
colormap jet;
for e = 1:num_elements
    nodes = elements1(e, corner_indices);
    patch(X(nodes), Y(nodes), stress_theta_nodal(nodes), 'FaceColor', 'interp', 'EdgeColor', 'none');
end
colorbar;
title('Hoop Stress Contour (\sigma_\theta)');
xlabel('X'); ylabel('Y');
axis equal; hold off;
