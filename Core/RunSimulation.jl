module RunSimulation

using ..GridSetup
using ..RFSolver
using ..BioheatSolver
using ..PlottingAndVisualization
using ..TimelapseCreation

export run_simulation

function run_simulation(grid_params, rf_params, bioheat_params, boundary_conditions)

    # Setup grid and solve RF problem
    grid = GridSetup.setup_grid(grid_params)
    V, dh, cellValues = solve_rf(grid, rf_params, grid_params, boundary_conditions)
    
    # Calculate E and Qel
    E = calculate_E(cellValues, dh, V)
    Qel, E_new = calculate_values(E, rf_params.sigma, grid_params)
    V_new = convert_V(V, rf_params.sigma, grid_params)



    #=
    # Initial temperature guess
    T_old = zeros(Float64, ndofs(dh))

    # Create timelapse or final visualization
    if create_timelapse
        TimelapseCreation.create_timelapse(grid, T_old, Qel, bioheat_params, grid_params)
    else
        # Solve bioheat equation without creating a timelapse
        T = BioheatSolver.solve_bioheat(grid, T_old, Qel, bioheat_params, grid_params)

        # Process and visualize the final temperature distribution
        PlottingAndVisualization.plot_slices(V, Qel, grid_params, reshape(V, grid_params.nx, grid_params.ny), reshape(Qel, grid_params.nx, grid_params.ny), T)
    end
    =#
    return grid, V, E, Qel, E_new, V_new
end

end  # module RunSimulation
