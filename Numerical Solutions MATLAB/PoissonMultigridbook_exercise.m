function PoissonSolver()

% Model parameters
radius = 6000;
density = 6000;
modelSize = 18000;
boundaryRadius = 8999;
numLevels = 4;
coarsestResolution = 7;
resolutionFactor = 2;
PoissonRelaxation = 1.5;
numSmoothingIterationsFinest = 5;
smoothingIterationFactor = 2;

% Create the grid
grid = CreateGrid(modelSize, boundaryRadius, numLevels, coarsestResolution, resolutionFactor);

% Initialize the potential
potential = zeros(grid.size);

% Do the smoothing iterations
for level = numLevels:-1:1
  numSmoothingIterations = numSmoothingIterationsFinest * smoothingIterationFactor ^ (level - 1);
  for i = 1:numSmoothingIterations
    potential = PoissonRelaxation * PoissonSolverStep(potential, grid) + (1 - PoissonRelaxation) * potential;
  end
end

% Write the potential to file
fid = fopen('potential.dat', 'w');
for i = 1:grid.size
  fprintf(fid, '%f\n', potential(i));
end
fclose(fid);

end

function grid = CreateGrid(modelSize, boundaryRadius, numLevels, coarsestResolution, resolutionFactor)

% Create the coarsest grid
coarsestGrid = Grid(modelSize, boundaryRadius, coarsestResolution);

% Create the finer grids
for level = 2:numLevels
  finerGrid = Grid(coarsestGrid.size, boundaryRadius, resolutionFactor);
  coarsestGrid = finerGrid;
end

% Combine the grids
grid = coarsestGrid;
for level = numLevels - 1:-1:1
  grid = CombineGrids(grid, finerGrid);
end

end

function potential = PoissonSolverStep(potential, grid)

% Compute the right-hand side of the Poisson equation
rhs = -density * grid.volume;

% Solve the Poisson equation
potential = PoissonSolverJacobi(rhs, grid);

end

function potential = PoissonSolverJacobi(rhs, grid)

% Initialize the potential
potential = zeros(grid.size);

% Iterate until convergence
for i = 1:100
  newPotential = PoissonRelaxation * rhs + (1 - PoissonRelaxation) * potential;
  if norm(newPotential - potential) < 1e-10
    break;
  end
  potential = newPotential;
end

end

function WritePotential(potential, grid)

% Write the potential to file
fid = fopen('potential.dat', 'w');
for i = 1:grid.size
  fprintf(fid, '%f\n', potential(i));
end
fclose(fid);

end

