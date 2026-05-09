# TimelapseCreation.jl

module TimelapseCreation
#=
using Plots
using FFMPEG

function create_timelapse(grid, T_old, Qel, bioheat_params, grid_params)
    nframes = bioheat_params.num_steps
    anim = @animate for step in 1:nframes
        # Update the temperature T for the current time step
        T = BioheatSolver.solve_bioheat(grid, T_old, Qel, bioheat_params, grid_params)
        
        # Update slice views for T
        img_t_x[1] = reshape(T[:, slice_y, :], grid_params.nx, grid_params.nz)
        img_t_y[1] = reshape(T[slice_x, :, :], grid_params.ny, grid_params.nz)
        img_t_z[1] = reshape(T[:, :, slice_z], grid_params.nx, grid_params.ny)

        # Plotting the temperature distribution as a heatmap
        heatmap(reshape(T, grid_params.nx, grid_params.ny), title="Time step: $step", color=:thermal)
        
        # Update T_old for the next iteration
        T_old = T
    end

    # Save the animation as a gif
    gif(anim, "heat_simulation_timelapse.gif", fps=10)
end
=#

end
