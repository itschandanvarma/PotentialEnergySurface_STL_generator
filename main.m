%% Potential Energy Surface STL Generator
% Author: Chandan Varma Tamada
clc; clear; close all;

% --- Input Parameters ---
% LEPS Parameters (energy in kcal/mol, distances in Angstroms):
% DFH: Dissociation energy of HF molecule (kcal/mol)
% DHH: Dissociation energy of H2 molecule (kcal/mol)
% betaFH: Morse potential exponent for HF bond (1/Angstrom), controls steepness of potential
% betaHH: Morse potential exponent for H2 bond (1/Angstrom), controls steepness of potential
% r0FH: Equilibrium bond length of HF (Angstrom)
% r0HH: Equilibrium bond length of H2 (Angstrom)
% SFH: Sato parameter for HF bond, adjusts interaction strength in LEPS model
% SHH: Sato parameter for H2 bond, adjusts interaction strength in LEPS model
DFH = 591.1;      % kcal/mol
DHH = 458.2;      % kcal/mol
betaFH = 2.2189;  % 1/Angstrom
betaHH = 1.9420;  % 1/Angstrom
r0FH = 0.917;     % Angstrom
r0HH = 0.7419;    % Angstrom
SFH = 0.167;
SHH = 0.106;

% Grid and Scaling Parameters:
% x_min_pes, x_max_pes: Range of r_FH (F-H distance) in Angstroms for PES grid
% y_min_pes, y_max_pes: Range of r_HH (H-H distance) in Angstroms for PES grid
% resolution: Number of grid points along each axis
% x_scale_factor, y_scale_factor: Scaling for STL model (mm/Angstrom)
% z_scale_factor: Scaling for potential in STL model (mm/(kcal/mol))
% z_limit_orig_units: Maximum energy (kcal/mol) to clip PES values
% floor_thickness_mm: Thickness of STL model base (mm)
% output_file_name: Name of output STL file
x_min_pes = 0.4;         % Angstroms
x_max_pes = 4.0;         % Angstroms
y_min_pes = 0.3;         % Angstroms
y_max_pes = 4.0;         % Angstroms
resolution = 600;        % Grid points
x_scale_factor = 40.0;   % mm / Angstrom
y_scale_factor = 40.0;   % mm / Angstrom
z_scale_factor = 0.05;   % mm / (kcal/mol)
z_limit_orig_units = 100; % kcal/mol
floor_thickness_mm = 2.0; % mm
output_file_name = 'PES.stl';

% --- Generate Grid ---
x_vals = linspace(x_min_pes, x_max_pes, resolution);
y_vals = linspace(y_min_pes, y_max_pes, resolution);
[X_orig, Y_orig] = meshgrid(x_vals, y_vals);

% --- Calculate LEPS PES ---
rac = X_orig + Y_orig;
den_ab = 1 + SFH;
den_ac = 1 + SFH;
den_bc = 1 + SHH;

% Q terms (vectorized)
Qab = 0.5 * DFH * (1.5 * exp(-2 * betaFH * (X_orig - r0FH)) - exp(-betaFH * (X_orig - r0FH))) / den_ab;
Qbc = 0.5 * DHH * (1.5 * exp(-2 * betaHH * (Y_orig - r0HH)) - exp(-betaHH * (Y_orig - r0HH))) / den_bc;
Qac = 0.5 * DFH * (1.5 * exp(-2 * betaFH * (rac - r0FH)) - exp(-betaFH * (rac - r0FH))) / den_ac;

% J terms (vectorized)
Jab = 0.5 * DFH * (0.5 * exp(-2 * betaFH * (X_orig - r0FH)) - 3 * exp(-betaFH * (X_orig - r0FH))) / den_ab;
Jbc = 0.5 * DHH * (0.5 * exp(-2 * betaHH * (Y_orig - r0HH)) - 3 * exp(-betaHH * (Y_orig - r0HH))) / den_bc;
Jac = 0.5 * DFH * (0.5 * exp(-2 * betaFH * (rac - r0FH)) - 3 * exp(-betaFH * (rac - r0FH))) / den_ac;

% LEPS potential
Z_surface = Qab + Qbc + Qac - sqrt(0.5 * ((Jab - Jbc).^2 + (Jbc - Jac).^2 + (Jac - Jab).^2));

% Clip and offset Z values
Z_surface = min(Z_surface, z_limit_orig_units); % Clip to max energy
Z_surface = Z_surface - min(Z_surface(:));      % Offset to zero

% --- 3D PES Plot ---
figure;
surf(X_orig, Y_orig, Z_surface, 'EdgeColor', 'none');
colormap('jet');
colorbar;
xlabel('r_{FH} (Å)');
ylabel('r_{HH} (Å)');
zlabel('Potential Energy');
title('LEPS Potential Energy Surface (kcal/mol)');
view(3);
shading interp;

% --- Generate STL File ---
% Scale coordinates
X_scaled = X_orig * x_scale_factor;
Y_scaled = Y_orig * y_scale_factor;
Z_scaled = Z_surface * z_scale_factor;
floor_z = -floor_thickness_mm;

% --- Create Vertices ---
n = resolution;
n_vertices = 2 * n * n; % Top surface + floor
vertices = zeros(n_vertices, 3);

% Top surface vertices
idx = 1:n^2;
vertices(idx,:) = [X_scaled(:), Y_scaled(:), Z_scaled(:)];

% Floor vertices
idx_floor = n^2+1:2*n^2;
vertices(idx_floor,:) = [X_scaled(:), Y_scaled(:), floor_z * ones(n^2,1)];

% --- Create Faces ---
% Preallocate faces array
n_faces = 2 * (n-1)^2 + 4 * 2 * (n-1); % Top + floor + 4 sides
faces = zeros(n_faces, 3); % Use double for triangulation compatibility

% Top surface triangles (normals up)
idx = 1:(n-1);
[J, I] = meshgrid(idx, idx);
v0 = (I-1)*n + J;
v1 = I*n + J;
v2 = (I-1)*n + (J+1);
v3 = I*n + (J+1);
face_idx = 1:2*(n-1)^2;
faces(face_idx,:) = [v0(:), v1(:), v2(:); v1(:), v3(:), v2(:)];

% Floor triangles (normals down)
f0 = n^2 + v0;
f1 = n^2 + v1;
f2 = n^2 + v2;
f3 = n^2 + v3;
face_idx = 2*(n-1)^2+1:4*(n-1)^2;
faces(face_idx,:) = [f0(:), f3(:), f1(:); f0(:), f2(:), f3(:)];

% Side walls
face_idx = 4*(n-1)^2+1;
% Left side (x_min, normal X-neg)
v_top = 1:n:(n-1)*n+1;
v_bot = n^2 + v_top;
faces(face_idx:face_idx+2*(n-1)-1,:) = [v_bot(1:end-1)', v_bot(2:end)', v_top(1:end-1)'; ...
                                         v_bot(2:end)', v_top(2:end)', v_top(1:end-1)'];
face_idx = face_idx + 2*(n-1);
% Right side (x_max, normal X-pos)
v_top = n:n:n^2;
v_bot = n^2 + v_top;
faces(face_idx:face_idx+2*(n-1)-1,:) = [v_top(1:end-1)', v_top(2:end)', v_bot(1:end-1)'; ...
                                         v_top(2:end)', v_bot(2:end)', v_bot(1:end-1)'];
face_idx = face_idx + 2*(n-1);
% Front side (y_min, normal Y-neg)
v_top = 1:n;
v_bot = n^2 + v_top;
faces(face_idx:face_idx+2*(n-1)-1,:) = [v_top(1:end-1)', v_bot(1:end-1)', v_bot(2:end)'; ...
                                         v_top(1:end-1)', v_bot(2:end)', v_top(2:end)'];
face_idx = face_idx + 2*(n-1);
% Back side (y_max, normal Y-pos)
v_top = (n-1)*n+1:n^2;
v_bot = n^2 + v_top;
faces(face_idx:face_idx+2*(n-1)-1,:) = [v_bot(1:end-1)', v_top(2:end)', v_top(1:end-1)'; ...
                                         v_bot(1:end-1)', v_bot(2:end)', v_top(2:end)'];

% --- Write STL File ---
TR = triangulation(double(faces), vertices); % Ensure faces is double
stlwrite(TR, output_file_name, 'binary');

% --- Visualize STL with Jet Colormap (No Edges) ---
figure;
trimesh(TR, 'FaceVertexCData', vertices(:,3), 'FaceColor', 'interp');
colormap('jet');
xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Z (mm)');
title('3D Visualization of LEPS PES STL Model');
axis equal;
grid on;
view(3);