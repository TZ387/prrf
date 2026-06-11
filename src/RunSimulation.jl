module RunSimulation

using ..GridSetup
using ..RFSolver
using ..BioheatSolver
using ..TimelapseCreation
using ..PlottingAndVisualization

export run_simulation, run_heat_simulation

function run_simulation(grid_params, rf_params, boundary_conditions)

    n_cells = grid_params.nx * grid_params.ny * grid_params.nz
    @info "Setting up grid ($(grid_params.nx)×$(grid_params.ny)×$(grid_params.nz), $n_cells cells)..."
    grid = GridSetup.setup_grid(grid_params)

    @info "Solving RF problem..."
    V_dof, dh, cellValues = solve_rf(grid, rf_params, grid_params, boundary_conditions)

    @info "Computing electric field and thermal dissipation..."
    Qel, E_mag, E_vec = calculate_fields(cellValues, dh, V_dof, rf_params.sigma, grid_params)
    V = convert_V(V_dof, grid_params)

    return grid, V_dof, Qel, E_mag, E_vec, V
end

function run_heat_simulation(Qel, grid_params, heat_params;
                        create_timelapse::Bool = false)

    T_final = nothing
    if create_timelapse
        @info "Running heat simulation with live plot..."
        T_final = TimelapseCreation.run_heat_timelapse(Qel, heat_params, grid_params)
    else
        @info "Running heat simulation (no live plot)..."
        Qzero = zeros(Float64, grid_params.nx, grid_params.ny, grid_params.nz)

        # Heating phase
        T = fill(heat_params.T_initial, grid_params.nx, grid_params.ny, grid_params.nz)
        T = BioheatSolver.solve_heat_phase(T, Qel, heat_params, grid_params,
                                            heat_params.t_on)
        # Cooling phase
        T_final = BioheatSolver.solve_heat_phase(T, Qzero, heat_params, grid_params,
                                                    heat_params.t_off)
    end

    return T_final
end

end  # module RunSimulation