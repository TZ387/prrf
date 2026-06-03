module RunSimulation

using ..GridSetup
using ..RFSolver
# using ..BioheatSolver
using ..PlottingAndVisualization
# using ..TimelapseCreation

export run_simulation

function run_simulation(grid_params, rf_params, bioheat_params, boundary_conditions)

    # Setup grid and solve RF problem
    grid = GridSetup.setup_grid(grid_params)
    V, dh, cellValues = solve_rf(grid, rf_params, grid_params, boundary_conditions)
    
    # Calculate E and Qel
    Qel, E_new = calculate_fields(cellValues, dh, V, rf_params.sigma, grid_params)
    V_new = convert_V(V, grid_params)

    return grid, V, Qel, E_new, V_new
end

end  # module RunSimulation
