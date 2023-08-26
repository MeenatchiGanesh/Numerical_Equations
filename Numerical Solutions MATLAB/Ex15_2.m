% Clear variables and close figures
clear all;
close all;

% Model parameters
radius = 6000; % Radius of the planetary body
density_planet = 6000; % Density of the planetary body
density_medium = 0; % Density of the surrounding medium
x_size = 18000; % Model size in x-direction
y_size = 18000; % Model size in y-direction
z_size = 18000; % Model size in z-direction
boundary_radius = 0.999 * x_size / 2; % Radius for boundary surface
level_num = 4; % Number of resolution levels
coarsest_resolution = 7; % Resolution on the coarsest grid
resolution_increase = 2; % Factor of increase in resolution between levels

% Multigrid parameters
relaxation_coefficient = 1.5; % Relaxation coefficient for Gauss-Seidel iterations

% Calculate number of grid points on each level
grid_points = zeros(level_num, 1);
grid_points(level_num) = coarsest_resolution;
for i = (level_num - 1):-1:1
    grid_points(i) = grid_points(i+1) * resolution_increase;
end

% Calculate grid steps on each level
grid_steps = x_size ./ (grid_points - 1);

% Initialize gravity potential phi on the finest grid
phi = zeros(grid_points(1), grid_points(1), grid_points(1));

% Main multigrid cycle
for level = 1:level_num
    
    % Set the number of smoothing iterations based on the current level
    num_iterations = 5 * (2^(level-1));
    
    % Perform smoothing iterations on the current level
    for iter = 1:num_iterations
        phi = Poisson3D_smoother(relaxation_coefficient, grid_points(level), grid_steps(level), phi, density_planet, density_medium, radius, boundary_radius);
    end
    
    % Restriction operation: Interpolate the gravity potential to the next coarser level
    if level < level_num
        phi_coarser = Poisson3D_restriction(phi, grid_points(level), grid_points(level+1));
        
        % Set the coarser level as the initial guess for the next iteration
        phi = phi_coarser;
    end
end

% Display the final gravity potential on the finest grid
% Clear variables and close figures
clear all;
close all;

% Model parameters
radius = 6000; % Radius of the planetary body
density_planet = 6000; % Density of the planetary body
density_medium = 0; % Density of the surrounding medium
x_size = 18000; % Model size in x-direction
y_size = 18000; % Model size in y-direction
z_size = 18000; % Model size in z-direction
boundary_radius = 0.999 * x_size / 2; % Radius for boundary surface
level_num = 4; % Number of resolution levels
coarsest_resolution = 7; % Resolution on the coarsest grid
resolution_increase = 2; % Factor of increase in resolution between levels

% Multigrid parameters
relaxation_coefficient = 1.5; % Relaxation coefficient for Gauss-Seidel iterations

% Calculate number of grid points on each level
grid_points = zeros(level_num, 1);
grid_points(level_num) = coarsest_resolution;
for i = (level_num - 1):-1:1
    grid_points(i) = grid_points(i+1) * resolution_increase;
end

% Calculate grid steps on each level
grid_steps = x_size ./ (grid_points - 1);

% Initialize gravity potential phi on the finest grid
phi = zeros(grid_points(1), grid_points(1), grid_points(1));

% Main multigrid cycle
for level = 1:level_num
    
    % Set the number of smoothing iterations based on the current level
    num_iterations = 5 * (2^(level-1));
    
    % Perform smoothing iterations on the current level
    for iter = 1:num_iterations
        phi = Poisson3D_smoother(relaxation_coefficient, grid_points(level), grid_steps(level), phi, density_planet, density_medium, radius, boundary_radius);
    end
    
    % Restriction operation: Interpolate the gravity potential to the next coarser level
    if level < level_num
        phi_coarser = Poisson3D_restriction(phi, grid_points(level), grid_points(level+1));
        
        % Set the coarser level as the initial guess for the next iteration
        phi = phi_coarser;
    end
end

% Display the final gravity potential 


slice(phi, [], [], (grid_points(1)-1)/2);
colormap jet;
shading interp;
colorbar;
xlabel('x');
ylabel('y');
zlabel('z');
title('Gravity Potential (Phi)');

% External function: Poisson3D_smoother
function phi_new = Poisson3D_smoother(relaxation_coefficient, grid_points, grid_step, phi, density_planet, density_medium, radius, boundary_radius)
    phi_new = phi; % Initialize the new gravity potential
    
    % Iterate over all grid points
    for i = 2:(grid_points-1)
        for j = 2:(grid_points-1)
            for k = 2:(grid_points-1)
                % Check if the grid point is inside the planetary body
                dist_from_center = grid_step * norm([i, j, k] - (grid_points+1)/2);
                if dist_from_center < radius
                    phi_new(i, j, k) = density_planet / density_medium * (phi(i+1, j, k) + phi(i-1, j, k) + phi(i, j+1, k) + phi(i, j-1, k) + phi(i, j, k+1) + phi(i, j, k-1)) / (6 * (grid_step^2));
                else
                    % Ghost node approach along the spherical boundary
                    if dist_from_center < boundary_radius
                        phi_new(i, j, k) = 0; % Set the ghost node potential to zero
                    else
                        % Compute the potential using the Laplace equation
                        phi_new(i, j, k) = (phi(i+1, j, k) + phi(i-1, j, k) + phi(i, j+1, k) + phi(i, j-1, k) + phi(i, j, k+1) + phi(i, j, k-1)) / 6;
                    end
                end
            end
        end
    end
end

% External function: Poisson3D_restriction
function phi_coarse = Poisson3D_restriction(phi_fine, grid_points_fine, grid_points_coarse)
    phi_coarse = zeros(grid_points_coarse, grid_points_coarse, grid_points_coarse);
    
    % Interpolate the gravity potential from the fine grid to the coarse grid
    for i = 1:grid_points_coarse
        for j = 1:grid_points_coarse
            for k = 1:grid_points_coarse
                i_fine = 2 * i - 1;
                j_fine = 2 * j - 1;
                k_fine = 2 * k - 1;
                
                phi_coarse(i, j, k) = (phi_fine(i_fine, j_fine, k_fine) + phi_fine(i_fine+1, j_fine, k_fine) + phi_fine(i_fine, j_fine+1, k_fine) + phi_fine(i_fine+1, j_fine+1, k_fine) ...
                    + phi_fine(i_fine, j_fine, k_fine+1) + phi_fine(i_fine+1, j_fine, k_fine+1) + phi_fine(i_fine, j_fine+1, k_fine+1) + phi_fine(i_fine+1, j_fine+1, k_fine+1)) / 8;
            end
        end
    end
end
