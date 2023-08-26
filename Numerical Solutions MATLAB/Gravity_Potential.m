radius_planet = 6000e3;  % Radius of the planet in meters
density_planet = 6000;   % Density of the planet in kg/m^3
model_size = 18000e3;    % Model size in meters
radius_boundary = 8999e3; % Radius for gravity potential boundary in meters

num_levels = 4;          % Number of resolution levels
coarsest_resolution = 7; % Resolution on the coarsest grid (last level)
resolution_increase = 2; % Factor of increase in resolution between the levels

relaxation_coefficient = 1.5; % Relaxation coefficient for Gauss-Seidel iterations
num_smoothing_iterations = 5; % Number of smoothing iterations on the finest level
iteration_increase = 2;       % Factor of increase in the number of iterations with level coarsening

grid_sizes = zeros(num_levels, 1);
grid_sizes(num_levels) = coarsest_resolution;
for level = num_levels-1:-1:1
    grid_sizes(level) = grid_sizes(level+1) * resolution_increase;
end

grid_spacing = model_size ./ grid_sizes;

ghost_node_boundary = 2;

gravity_potential = zeros(grid_sizes(1), grid_sizes(1), grid_sizes(1));
for level = num_levels:-1:1
    gravity_potential = computeGravityPotential(level, grid_sizes, grid_spacing, gravity_potential);
end

[X, Y, Z] = meshgrid(1:grid_sizes(1), 1:grid_sizes(1), 1:grid_sizes(1));
figure;
isosurface(X, Y, Z, gravity_potential);
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Gravity Potential');

function potential = computeGravityPotential(level, grid_sizes, grid_spacing, gravity_potential)
    grid_size = grid_sizes(level);
    spacing = grid_spacing(level);
    relaxation = relaxation_coefficient;
    num_iterations = num_smoothing_iterations * (iteration_increase^(num_levels-level));

    % Perform the computations for the current level
    potential = zeros(grid_size, grid_size, grid_size);

    % Set up the grid points
    x = linspace(-model_size/2, model_size/2, grid_size);
    y = linspace(-model_size/2, model_size/2, grid_size);
    z = linspace(-model_size/2, model_size/2, grid_size);
    [X, Y, Z] = meshgrid(x, y, z);
    
    % Compute the gravity potential for Mars-like planetary body
    r = sqrt(X.^2 + Y.^2 + Z.^2);
    r_planet = radius_planet * ones(size(r));
    r_boundary = radius_boundary * ones(size(r));
    rho = density_planet * ones(size(r));
    rho(r > r_boundary) = 0;
    
    for iteration = 1:num_iterations
        for i = 2:grid_size-1
            for j = 2:grid_size-1
                for k = 2:grid_size-1
                    if r(i, j, k) <= r_planet
                        potential(i, j, k) = (1/6) * (potential(i-1, j, k) + potential(i+1, j, k) ...
                            + potential(i, j-1, k) + potential(i, j+1, k) + potential(i, j, k-1) ...
                            + potential(i, j, k+1) - (rho(i, j, k)/epsilon_0) * spacing^2);
                    end
                end
            end
        end
    end

    % Apply the ghost node approach along the spherical boundary
    potential(1, :, :) = -potential(2, :, :);
    potential(grid_size, :, :) = -potential(grid_size-1, :, :);
    potential(:, 1, :) = -potential(:, 2, :);
    potential(:, grid_size, :) = -potential(:, grid_size-1, :);
    potential(:, :, 1) = -potential(:, :, 2);
    potential(:, :, grid_size) = -potential(:, :, grid_size-1);

    potential = potential + potential(1, 1, 1);

    % Return the computed gravity potential at the current level
    potential = potential .* (1e-9);  % Scale the potential for visualization
end
