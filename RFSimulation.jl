# RFSimulation.jl - Main module that includes all submodules

module RFSimulation

# Include all the core modules
include("Core/Config.jl")
include("Core/GridSetup.jl") 
include("Core/RFSolver.jl")
# include("Core/BioheatSolver.jl")
include("Core/PlottingAndVisualization.jl")
# include("Core/TimelapseCreation.jl")
include("Core/RunSimulation.jl")

# Re-export the submodules so they can be accessed directly
using .Config
using .GridSetup
using .RFSolver
# using .BioheatSolver
using .PlottingAndVisualization
# using .TimelapseCreation
using .RunSimulation

# Optionally, you can also re-export specific functions/types for convenience
export GridParams, RFParams, BioheatParams, setup_material_properties, create_coordinate_grids
export setup_grid
export solve_rf, calculate_E, calculate_values, convert_V
export plot_slices, plot_graphs
export run_simulation

end # module RFSimulation