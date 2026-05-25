module prrf

include("Config.jl")
include("GridSetup.jl")
include("RFSolver.jl")
include("PlottingAndVisualization.jl")
include("RunSimulation.jl")

using .Config
using .GridSetup
using .RFSolver
using .PlottingAndVisualization
using .RunSimulation

export GridParams, RFParams, BioheatParams, setup_material_properties, create_coordinate_grids
export setup_grid
export solve_rf, calculate_E, calculate_values, convert_V
export plot_slices, plot_graphs
export run_simulation

end