% Clearing all variables and arrays
clear all;
% Clearing figures
clf;

% Model parameters
% Radius of the planet
r = 6000000.0; % 6000 km
% Density of the planet
rhoplanet = 6000.0; % 6000 kg/m^3
% Gravity constant
G = 6.67e-11;
% Model size, m
xsize = 18000000.0; % 18000 km
ysize = 18000000.0; % 18000 km
zsize = 18000000.0; % 18000 km
% Radius for gravity potential boundary
gradius = 8999000.0; % 8999 km

% Multigrid parameters
% Number of resolution levels
levelnum = 4;

% Iteration Parameters
% Total number of iteration cycles
inum = 30;
% Relaxation coefficient for Gauss-Seidel iterations
krelax = 1.5;
% Number of smoothing iterations
iternum(1) = 5;
iternum(2) = 5 * 2;
iternum(3) = 5 * 4;
iternum(4) = 5 * 8;
iternum(5) = 5 * 16;
iternum(6) = 5 * 32;
iternum(7) = 5 * 64;
iternum(8) = 5 * 128;

% Defining resolutions on all levels
xnum(1) = 49;
ynum(1) = 49;
znum(1) = 49;
xnum(2) = 25;
ynum(2) = 25;
znum(2) = 25;
xnum(3) = 13;
ynum(3) = 13;
znum(3) = 13;
xnum(4) = 7;
ynum(4) = 7;
znum(4) = 7;

% Defining gridsteps on all levels
for n = 1:levelnum
    xstp(n) = xsize / (xnum(n) - 1);
    ystp(n) = ysize / (ynum(n) - 1);
    zstp(n) = zsize / (znum(n) - 1);
end

% FINEST (PRINCIPAL) GRID
% Defining density structure rho()
% Defining initial guesses for gravity potential phi()
% Computing right part of the Poisson equation R()
% Grid points cycle
for i = 1:ynum(1)
    for j = 1:xnum(1)
        for k = 1:znum(1)

            % Gravity potential
            phi1(i, j, k) = 0;
            % Check distance of (i,j) node from the grid center
            dx = (j - 1) * xstp(1) - xsize / 2;
            dy = (i - 1) * ystp(1) - ysize / 2;
            dz = (k - 1) * zstp(1) - zsize / 2;
            dist = sqrt(dx * dx + dy * dy + dz * dz);
            % Density
            rho(i, j, k) = 0;
            if (dist < r)
                rho(i, j, k) = rhoplanet;
            end
            % Right part of Poisson equation
            R1(i, j, k) = 4.0 * pi * G * rhoplanet;
        end
    end
end

% Defining boundary condition nodes for all grids
for n = levelnum:-1:1
    % Grid points cycle
    for i = 1:ynum(n)
        for j = 1:xnum(n)
            for k = 1:xnum(n)
                % Check distance of (i,j) node from the grid center
                dx = (j - 1) * xstp(n) - xsize / 2;
                dy = (i - 1) * ystp(n) - ysize / 2;
                dz = (k - 1) * zstp(n) - zsize / 2;
                dist = sqrt(dx * dx + dy * dy + dz * dz);
                % Action depends on the level of resolution
                switch n
                    case 1
                        bon1(i, j, k) = 0;
                        if (dist < gradius)
                            bon1(i, j, k) = 1;
                        end
                    case 2
                        bon2(i, j, k) = 0;
                        if (dist < gradius)
                            bon2(i, j, k) = 1;
                        end
                    case 3
                        bon3(i, j, k) = 0;
                        if (dist < gradius)
                            bon3(i, j, k) = 1;
                        end
                    case 4
                        bon4(i, j, k) = 0;
                        if (dist < gradius)
                            bon4(i, j, k) = 1;
                        end
                end
            end
        end
    end
end

% Figures counter
fignum = 1;

% Main Multigrid cycle
for niter = 1:inum

    % Smoothing+restriction cycle
    for n = 1:levelnum
        % Action depends on the level of resolution
        switch n

            % Level 1 (principal grid)
            case 1
                % Smoothing operation: solving of Poisson equation on nodes
                % and computing residuals
                [phi1, residual1] = Poisson3D_smoother_planet(iternum(n), krelax, xnum(n), ynum(n), znum(n), xstp(n), ystp(n), zstp(n), R1, phi1, bon1, gradius);
                % Restriction operation:
                % Interpolating residuals to coarser level (k+1)
                % to produce right parts for this coarser level
                if (levelnum > n)
                    [R2] = Poisson3D_restriction_planet(n, xnum, ynum, znum, xstp, ystp, zstp, residual1, bon1, bon2);
                end

            % Level 2
            case 2
                % Making initial approximation for phi corrections (zeros)
                phi2 = zeros(ynum(n), xnum(n), znum(n));
                % Smoothing operation:
                [phi2, residual2] = Poisson3D_smoother_planet(iternum(n), krelax, xnum(n), ynum(n), znum(n), xstp(n), ystp(n), zstp(n), R2, phi2, bon2, gradius);
                % Restriction operation:
                if (levelnum > n)
                    [R3] = Poisson3D_restriction_planet(n, xnum, ynum, znum, xstp, ystp, zstp, residual2, bon2, bon3);
                end

            % Level 3
            case 3
                % Making initial approximation for phi corrections (zeros)
                phi3 = zeros(ynum(n), xnum(n), znum(n));
                % Smoothing operation:
                [phi3, residual3] = Poisson3D_smoother_planet(iternum(n), krelax, xnum(n), ynum(n), znum(n), xstp(n), ystp(n), zstp(n), R3, phi3, bon3, gradius);
                % Restriction operation:
                if (levelnum > n)
                    [R4] = Poisson3D_restriction_planet(n, xnum, ynum, znum, xstp, ystp, zstp, residual3, bon3, bon4);
                end

            % Level 4
            case 4
                % Making initial approximation for phi corrections (zeros)
                phi4 = zeros(ynum(n), xnum(n), znum(n));
                % Smoothing operation:
                [phi4, residual4] = Poisson3D_smoother_planet(iternum(n), krelax, xnum(n), ynum(n), znum(n), xstp(n), ystp(n), zstp(n), R4, phi4, bon4, gradius);
        end
    end

    % Plotting Residuals for selected stages
    % Normalizing residual
    kfr = 4.0 * pi * G * rhoplanet;
    residual0 = residual1 / kfr;

    % Plotting Potential
    figure(1);
    subplot(2, 2, 1);
    surf(phi1(:, :, (znum(1) - 1) / 2));
    shading interp;
    light;
    lighting phong;
    axis tight;
    zlabel('Gravity potential');
    title(['Solution of Poisson equation, V-cycle = ', num2str(niter)]);

    % Plotting Residual
    subplot(2, 2, 2);
    surf(residual0(:, :, (znum(1) - 1) / 2));
    shading interp;
    light;
    lighting phong;
    axis tight;
    zlabel('Residual');
    title(['Solution of Poisson equation, V-cycle = ', num2str(niter)]);

    % Computing mean square residuals
    resphi00(niter) = sum(sum(sum(residual0 .^ 2))) / (ynum(1) * xnum(1) * znum(1));
    % Compute log of new mean square residuals
    resphi00(niter) = log10(sqrt(resphi00(niter)));

    % Plotting Mean residuals
    subplot(2, 2, 3);
    plot(resphi00, 'k');
    xlabel('V-cycles');
    ylabel(['log(Residuals) iteration ', num2str(niter)]);

    pause(0.5);
end
