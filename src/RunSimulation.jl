module RunSimulation

using ..GridSetup
using ..RFSolver
# using ..BioheatSolver
using ..PlottingAndVisualization
# using ..TimelapseCreation

export run_simulation

function run_simulation(grid_params, rf_params, bioheat_params, boundary_conditions)
    n_cells = grid_params.nx * grid_params.ny * grid_params.nz
    @info "Setting up grid ($(grid_params.nx)×$(grid_params.ny)×$(grid_params.nz), $n_cells cells)..."
    grid = GridSetup.setup_grid(grid_params)

    @info "Solving RF problem..."
    V_dof, dh, cellValues = solve_rf(grid, rf_params, grid_params, boundary_conditions)

    @info "Computing electric field and thermal dissipation..."
    Qel, E = calculate_fields(cellValues, dh, V_dof, rf_params.sigma, grid_params)
    V = convert_V(V_dof, grid_params)

    return grid, V_dof, Qel, E, V
end

end  # module RunSimulation
