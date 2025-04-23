    clc;
clearvars;
clf;
%% Parameters   
ri = input("Enter internal radius="); % Inner radius of the cylinder (m)
ro = input("Enter external radius="); % outer radius of the cylinder (m)
L = input("Enter lenght of the cylinder="); % lenght of the cylinder (m)
nr = input("enter number of divisions in r direction=");  % Number of elements along r direction
nz = input("enter numbe of divisions in z direction=");  % Number of elements along y

E = input("enter Youngs moduLus in MPa=");   % Young's modulus (Pa)
v = input("enter poissons ratio=");     % Poisson's ratio
P = input("enter internal pressue in MPa= ");    % Pressure applied on the right side (N/m²)

% Constitutive matrix
C = (E / ((1+v)*(1-2*v))) * [1-v, v, v, 0; 
                             v, 1-v, v, 0; 
                             v, v, 1-v, 0;
                             0, 0, 0, (1-2*v)/2];

%% Generate Node Coordinates
r_vals = linspace(ri, ro, nr + 1);  % Radial nodes from inner to outer radius
z_vals = linspace(0, L, nz + 1);    % Axial nodes from 0 to height H

node_coords = [];
for j = 1:length(z_vals)
    for i = 1:length(r_vals) 
        node_coords = [node_coords; r_vals(i), z_vals(j)];
    end
end

num_nodes = size(node_coords, 1);

%% Generate Elements (4-node Quadrilaterals)
elements = [];
for j = 1:nz  % Axial direction (Z)
    for i = 1:nr  % Radial direction (R)
        n1 = (j - 1) * (nr + 1) + i;
        n2 = n1 + 1;
        n3 = n2 + (nr + 1);
        n4 = n1 + (nr + 1);
        elements = [elements; n1, n2, n3, n4];
    end
end

num_elements = size(elements, 1);


%% Visualization of Mesh
figure;
hold on;
for elem = 1:num_elements
    nodes = elements(elem, :);
    coords = node_coords(nodes, :);
    fill(coords(:,1), coords(:,2), 'w', 'EdgeColor', 'k'); % Plot elements
end
axis equal;
title('Axisymmetric Mesh (r-z Plane)');
xlabel('Radial Coordinate (r)');
ylabel('Axial Coordinate (z)');
%% Gauss Quadrature Points and Weights (2x2)
gauss_points = [-1/sqrt(3), 1/sqrt(3)];
weights = [1, 1];
%% Define Symbolic Shape Functions
syms xi eta

% Shape functions for 4-node quadrilateral
N1 = (1 - xi) * (1 - eta) / 4;
N2 = (1 + xi) * (1 - eta) / 4;
N3 = (1 + xi) * (1 + eta) / 4;
N4 = (1 - xi) * (1 + eta) / 4;
Ni = [N1, N2, N3, N4];

% Compute derivatives symbolically
dN_dxi = [eta/4 - 1/4, 1/4 - eta/4, eta/4 + 1/4, - eta/4 - 1/4];
dN_deta = [xi/4 - 1/4, - xi/4 - 1/4, xi/4 + 1/4, 1/4 - xi/4]; 
%% Initialize Global Stiffness Matrix
K_global = zeros(2 * num_nodes, 2 * num_nodes);
F_global = zeros(2 * num_nodes, 1);
%% Compute Stiffness Matrix for Axisymmetric Analysis
for elem = 1:num_elements
    % Extract node indices and coordinates for the current element
    node_ids = elements(elem, :);
    coords = node_coords(node_ids, :);
    
    % Element stiffness matrix
    Ke = zeros(8, 8);
    
    % Gauss Integration
    for i = 1:2
        for j = 1:2
            xi_val = gauss_points(i);
            eta_val = gauss_points(j);
            
            % Evaluate shape function derivatives at Gauss points
            dN_dxi_val = double(subs(dN_dxi, [xi, eta], [xi_val, eta_val]));
            dN_deta_val = double(subs(dN_deta, [xi, eta], [xi_val, eta_val]));

            % Compute Jacobian Matrix
            J = [dN_dxi_val; dN_deta_val] * coords;
            detJ = det(J);
            Jinv = inv(J);
            
            % Compute B matrix for axisymmetric analysis
            B = zeros(4, 8);
            
            % Compute r at Gauss point correctly
            Ni_val = double(subs(Ni, [xi, eta], [xi_val, eta_val]));
            r_gauss = Ni_val * coords(:,1);
            
            for k = 1:4
                dN_nat = [dN_dxi_val(k); dN_deta_val(k)];
                dN_xy = Jinv * dN_nat;
                
                B(1, 2*k-1) = dN_xy(1);  % dN/dr
                B(2, 2*k)   = dN_xy(2);  % dN/dz
                B(3, 2*k-1) = Ni_val(k) / r_gauss;  % Hoop strain term (N/r)
                B(4, 2*k-1) = dN_xy(2);  % dN/dz
                B(4, 2*k)   = dN_xy(1);  % dN/dr
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
%% Apply Boundary Conditions (Fix Only Uz at z = 0 and z = L)
fixed_nodes_z = find(node_coords(:,2) == 0 | node_coords(:,2) == L); % Nodes at z = 0 (bottom) and z = L (top)
fixed_dofs_z = 2 * fixed_nodes_z; % Fix only Uz (longitudinal displacement)

fixed_dofs = unique(fixed_dofs_z(:)); % Ensure unique DOFs

% Apply boundary conditions to stiffness and force matrices
K_global(fixed_dofs, :) = 0; 
K_global(:, fixed_dofs) = 0; 
K_global(fixed_dofs, fixed_dofs) = eye(length(fixed_dofs)); 
F_global(fixed_dofs) = 0; % Ensure force vector consistency

%% Apply Internal Pressure Load at Inner Radius (r = ri)
inner_nodes = find(node_coords(:,1) == ri); % Find nodes at inner radius
inner_dofs = 2 * inner_nodes - 1; % Ur (radial displacement) DOFs

% Compute total force per unit length (Pressure × circumference × element length)
element_length = L / nz;  % Length of each element along z
total_force = P * 2 * pi * ri * element_length;  % Total force per element (axisymmetric)

% Distribute force across inner-edge nodes using shape functions
for i = 1:length(inner_nodes)
    if i == 1  % First node (only connected to one element)
        F_global(inner_dofs(i)) = total_force / 2;  % Apply half force for first node
    elseif i == length(inner_nodes)  % Last node (only connected to one element)
        F_global(inner_dofs(i)) = total_force / 2;  % Apply half force for last node
    else  % Internal nodes (connected to two elements)
        F_global(inner_dofs(i)) = total_force;  % Full force for internal nodes
    end
end

%% Solve for Displacements

% Ensure no singularities by modifying fixed DOFs
nonzero_dofs = setdiff(1:length(F_global), fixed_dofs); % Remove fixed DOFs

% Extract submatrices for non-fixed DOFs
K_reduced = K_global(nonzero_dofs, nonzero_dofs);
F_reduced = F_global(nonzero_dofs);

% Solve for reduced displacements
U_reduced = K_reduced \ F_reduced;

% Initialize full displacement vector
U = zeros(length(F_global), 1);
U(nonzero_dofs) = U_reduced; % Assign computed values


%% Extract Displacement Results
Ux = U(1:2:end); % Radial displacements (Ur)
Uz = U(2:2:end); % Axial displacements (Uz)
max(Ux)
min(Ux)
min(Uz)
%% Generate Grid Data
x = node_coords(:,1); % Radial coordinate (r)
y = node_coords(:,2); % Axial coordinate (z)

% Create interpolant for radial displacement (Ux)
Fx = scatteredInterpolant(x, y, Ux, 'natural', 'none'); 

% Generate Regular Grid for Contour Plot
[xq, yq] = meshgrid(linspace(min(x), max(x), 100), linspace(min(y), max(y), 100));
Uxq = Fx(xq, yq);  % Interpolated radial displacement

%% Plot Radial Displacement Contour (Ur)
figure;
contourf(xq, yq, Uxq, 30, 'LineColor', 'none'); % Smooth contour plot
colorbar;
title('Radial Displacement Contour (Ur)');
xlabel('Radial Coordinate (r)');
ylabel('Axial Coordinate (z)');
axis equal;

%% Compute Stress at Each Element
element_stress = zeros(num_elements, 4); % Stresses (σ_rr, σ_zz, σ_θθ, τ_rz)

for elem = 1:num_elements
    node_ids = elements(elem, :);
    coords = node_coords(node_ids, :);
    
    % Extract displacement values for this element (Fixed)
    dof_map = reshape([2*node_ids-1; 2*node_ids], [], 1); % Ensure column vector
    Ue = U(dof_map);
    
    % Initialize stress accumulator
    stress_at_gauss = zeros(4, 4); % (4 stress components) × (4 Gauss points)

    for i = 1:2
        for j = 1:2
            xi_val = gauss_points(i);
            eta_val = gauss_points(j);
            
            % Compute shape function derivatives
            dN_dxi_val = double(subs(dN_dxi, [xi, eta], [xi_val, eta_val]));
            dN_deta_val = double(subs(dN_deta, [xi, eta], [xi_val, eta_val]));
            
            % Compute Jacobian and inverse
            J = [dN_dxi_val; dN_deta_val] * coords;
            Jinv = inv(J);
            
            % Compute B-matrix
            B = zeros(4, 8);
            Ni_val = double(subs(Ni, [xi, eta], [xi_val, eta_val]));
            r_gauss = Ni_val * coords(:,1);
            
            for k = 1:4
                dN_nat = [dN_dxi_val(k); dN_deta_val(k)];
                dN_xy = Jinv * dN_nat;
                
                B(1, 2*k-1) = dN_xy(1);  % dN/dr
                B(2, 2*k)   = dN_xy(2);  % dN/dz
                B(3, 2*k-1) = Ni_val(k) / r_gauss;  % Hoop strain term
                B(4, 2*k-1) = dN_xy(2);  % dN/dz
                B(4, 2*k)   = dN_xy(1);  % dN/dr
            end
            
            % Compute stress at Gauss point
            stress_at_gauss(:, (i-1)*2 + j) = C * B * Ue;
        end
    end
    
    % Average stress over Gauss points
    element_stress(elem, :) = mean(stress_at_gauss, 2);
end

%% Convert Elemental Stresses to Nodal Stresses
nodal_stress_rr = zeros(num_nodes, 1);
nodal_stress_zz = zeros(num_nodes, 1);
node_count = zeros(num_nodes, 1);

for elem = 1:num_elements
    node_ids = elements(elem, :);
    
    % Accumulate stress values
    nodal_stress_rr(node_ids) = nodal_stress_rr(node_ids) + element_stress(elem, 1); % σ_rr
    nodal_stress_zz(node_ids) = nodal_stress_zz(node_ids) + element_stress(elem, 2); % σ_zz
    node_count(node_ids) = node_count(node_ids) + 1;
end

% Average over shared nodes
nodal_stress_rr = nodal_stress_rr ./ node_count;
minimum_radial_stress=min(nodal_stress_rr)
maximum_radial_stress=max(nodal_stress_rr)
nodal_stress_zz = nodal_stress_zz ./ node_count;
minimum_axial_stress=min(nodal_stress_zz)
maximum_axial_stress=max(nodal_stress_zz)
%% Plot Radial Stress (σ_rr)
figure;
patch('Faces', elements, 'Vertices', node_coords, ...
      'FaceVertexCData', nodal_stress_rr, 'EdgeColor', 'k', 'FaceColor', 'interp');
colorbar;
title('Radial Stress (\sigma_{rr}) Distribution');
xlabel('Radial Coordinate (r)');
ylabel('Axial Coordinate (z)');
axis equal;

%% Plot Axial Stress (σ_zz)
figure;
patch('Faces', elements, 'Vertices', node_coords, ...
      'FaceVertexCData', nodal_stress_zz, 'EdgeColor', 'k', 'FaceColor', 'interp');
colorbar;
title('Axial Stress (\sigma_{zz}) Distribution');
xlabel('Radial Coordinate (r)');
ylabel('Axial Coordinate (z)');
axis equal;

%% Compute Hoop Stress (σ_θθ) at Nodes
nodal_stress_theta = zeros(num_nodes, 1);

for elem = 1:num_elements
    node_ids = elements(elem, :);
    
    % Accumulate hoop stresses
    nodal_stress_theta(node_ids) = nodal_stress_theta(node_ids) + element_stress(elem, 3); % σ_θθ
end

% Average over shared nodes
nodal_stress_theta = nodal_stress_theta ./ node_count;

% Find min/max values
min_hoop_stress = min(nodal_stress_theta);
max_hoop_stress = max(nodal_stress_theta);
disp(['Minimum Hoop Stress (σ_θθ_min): ', num2str(min_hoop_stress)]);
disp(['Maximum Hoop Stress (σ_θθ_max): ', num2str(max_hoop_stress)]);

%% Plot Hoop Stress (σ_θθ) Distribution
figure;
patch('Faces', elements, 'Vertices', node_coords, ...
      'FaceVertexCData', nodal_stress_theta, 'EdgeColor', 'k', 'FaceColor', 'interp');
colorbar;
title('Hoop Stress (\sigma_{\theta\theta}) Distribution');
xlabel('Radial Coordinate (r)');
ylabel('Axial Coordinate (z)');
axis equal;