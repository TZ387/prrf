# Include the main module - this loads all submodules
using prrf

# Geometry definition function
function geometryDefinition(grid_params::GridParams)
    # Get coordinate grids from Config module (implementation details hidden)
    X, Y, Z = create_coordinate_grids(grid_params)
    # Initialize with default material (air/gel)
    M = ones(Int, grid_params.nx, grid_params.ny, grid_params.nz)
    
    # Apply conditions using absolute coordinates [m]
    M[Z .> 0.0005] .= 2    # Skin (0.5 mm depth)
    M[Z .> 0.002] .= 3     # Fat (2 mm depth)  
    M[Z .> 0.0055] .= 4    # Muscle (5.5 mm depth)
    
    return M
end

# Define materials and their properties
function define_material_properties()
    materials = Dict{Int, NamedTuple}()
    
    epsilon_0 = 8.854e-12  # Permittivity of free space [F/m]
    materials[1] = (name = "Gel", sigma = 0.3, epsilon_im = 0*epsilon_0, VHC = 1.0e6, k = 0.025)
    materials[2] = (name = "Skin", sigma = 0.22, epsilon_im = 0*epsilon_0, VHC = 4.18e6, k = 0.6)
    materials[3] = (name = "Fat", sigma = 0.025, epsilon_im = 0*epsilon_0, VHC = 1.0e6, k = 0.5)
    materials[4] = (name = "Muscle", sigma = 0.5, epsilon_im = 0*epsilon_0, VHC = 1.0e6, k = 0.5)
    
    return materials
end

function main()
    # Flags for use
    use_gpu = false              # Decide if you want to use GPU for calculations or not
    create_timelapse = false     # Decide if you want to create a timelapse of the heat simulation

    # Define grid parameters
    grid_params = GridParams(
        0.3,    # lx = Length of the simulation domain in x-direction [m]
        0.3,    # ly = Length of the simulation domain in y-direction [m]
        0.0105, # lz = Length of the simulation domain in z-direction [m]
        100,    # nx = Number of elements in the x-direction
        100,    # ny = Number of elements in the y-direction
        210     # nz = Number of elements in the z-direction
    )

    # Load materials and geometry
    materials = define_material_properties()
    material_indices = geometryDefinition(grid_params)

    # Setup material properties (delegated to Config.jl)
    sigma, epsilon_im, VHC, k = setup_material_properties(material_indices, materials, grid_params.nx, grid_params.ny, grid_params.nz)

    # Load RF parameters
    rf_params = RFParams(
        sigma,   # Electrical conductivity matrix [S/m]
        epsilon_im, # Permittivity matrix [F/m]
        1e6,     # ω = (Circular) Frequency of the RF signal [Hz]
    )

    # Load bioheat parameters
    bioheat_params = BioheatParams(
        1.0,     # Δt = Time step size [s]
        100,     # num_steps = Number of time steps to simulate
        VHC,       # Volumetric heat capacity matrix [J/(m³·K)]
        k,       # Thermal conductivity matrix [W/(m·K)]
        1050.0,  # rho_b = Blood density [kg/m³]
        3600.0,  # cb = Specific heat capacity of blood [J/(kg·K)]
        0.001,   # ωb = Blood perfusion rate [1/s]
        37.0,    # Ta = Arterial temperature [°C]
        0.0      # Qmet = Metabolic heat generation term [W/m³]
    )

    # Define boundary conditions as functions
    boundary_conditions = Dict(
        "bottom" => (x, t) -> 10.0,
        "top" => (x, t) -> 0.0
        # Add other boundaries here as needed
    )


    # Run the simulation using the RunSimulation module and plot graphs
    grid, V, E, Qel, E_new, V_new = run_simulation(grid_params, rf_params, bioheat_params, boundary_conditions);

    # Save
    save_simulation("Example1.h5", grid_params, material_indices;
        Qel=Qel, E_new=E_new, V_new=V_new)

    plot_graphs(material_indices, grid_params, Qel, E_new, V_new)

    return grid, V, E, Qel, E_new, V_new
end

# Call the main function and return values for potential inspection
grid, V, E, Qel, E_new, V_new = main();