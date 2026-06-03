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
    V_dof, dh, cellValues = solve_rf(grid, rf_params, grid_params, boundary_conditions)
    
    # Calculate E and Qel
    Qel, E = calculate_fields(cellValues, dh, V_dof, rf_params.sigma, grid_params)
    V = convert_V(V_dof, grid_params)

    return grid, V_dof, Qel, E, V
end

end  # module RunSimulation
