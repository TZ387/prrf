# Config.jl

module Config

# Export the structures and functions to make them available outside this module
export RFParams, BioheatParams, GridParams, setup_material_properties, create_coordinate_grids

# Parameters for RF problem
struct RFParams
    sigma::Array{Float64, 3}  # Electrical conductivity matrix [S/m]
    epsilon_im::Array{Float64, 3}  # Permittivity matrix [F/m]
    ω::Float64  # (Circular) Frequency of the RF signal [Hz]
end

# Parameters for bioheat problem
struct BioheatParams
    Δt::Float64  # Time step size [s]
    num_steps::Int  # Number of time steps to simulate
    VHC::Array{Float64, 3}  # Volumetric heat capacity matrix [J/(m³·K)] (VHC = rho * c)
    k::Array{Float64, 3}  # Thermal conductivity matrix [W/(m·K)]
    rho_b::Float64  # Blood density [kg/m³]
    cb::Float64  # Specific heat capacity of blood [J/(kg·K)]
    ωb::Float64  # Blood perfusion rate [1/s]
    Ta::Float64  # Arterial temperature [°C]
    Qmet::Float64  # Metabolic heat generation term [W/m³]
end

# Parameters for grid setup
struct GridParams
    lx::Float64  # Length of the simulation domain in x-direction [m]
    ly::Float64  # Length of the simulation domain in y-direction [m]
    lz::Float64  # Length of the simulation domain in z-direction [m]
    nx::Int  # Number of elements in the x-direction
    ny::Int  # Number of elements in the y-direction
    nz::Int  # Number of elements in the z-direction
end

# Function to create 3D coordinate grids from grid parameters
function create_coordinate_grids(grid_params::GridParams)
    # Create 3D coordinate grids - contains actual coordinates in meters
    # Use (index - 0.5) to get cell center coordinates
    X = repeat(reshape((1:grid_params.nx) .- 0.5, grid_params.nx, 1, 1) .* (grid_params.lx / grid_params.nx), 1, grid_params.ny, grid_params.nz)
    Y = repeat(reshape((1:grid_params.ny) .- 0.5, 1, grid_params.ny, 1) .* (grid_params.ly / grid_params.ny), grid_params.nx, 1, grid_params.nz)
    Z = repeat(reshape((1:grid_params.nz) .- 0.5, 1, 1, grid_params.nz) .* (grid_params.lz / grid_params.nz), grid_params.nx, grid_params.ny, 1)
    
    return X, Y, Z
end

# Function to setup material properties matrices based on the geometry
function setup_material_properties(material_indices, materials, nx, ny, nz)
    sigma = zeros(Float64, nx, ny, nz)
    epsilon_im = zeros(Float64, nx, ny, nz)
    VHC = zeros(Float64, nx, ny, nz)
    k = zeros(Float64, nx, ny, nz)
    
    for i in eachindex(material_indices)
        mat = materials[material_indices[i]]
        sigma[i]      = mat[:sigma]
        epsilon_im[i] = mat[:epsilon_im]
        VHC[i]        = mat[:VHC]
        k[i]          = mat[:k]
    end
    
    return sigma, epsilon_im, VHC, k
end

end # module Config