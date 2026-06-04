module PlottingAndVisualization

using GLMakie
using ..Config

export plot_slices, plot_graphs

function plot_slices(V, grid_params; fig = nothing, num_ticks = 5, use_log_scale = false, title = "")
    if fig === nothing
        fig = Figure(size = (700, 700))
    end

    ax = LScene(fig[1, 1], show_axis=false)
    maximum_length = max(grid_params.lx, grid_params.ly, grid_params.lz)
    scale!(ax.scene, maximum_length/grid_params.lx, maximum_length/grid_params.ly, maximum_length/grid_params.lz)

    # Set the title for the figure
    if !isempty(title)
        Label(fig[0, 1], title, fontsize=18, tellwidth=false)  # Add title at the top
    end

    x = LinRange(-grid_params.lx/2, grid_params.lx/2, grid_params.nx)
    y = LinRange(-grid_params.ly/2, grid_params.ly/2, grid_params.ny)
    z = LinRange(grid_params.lz, 0, grid_params.nz)

    if use_log_scale
        # Apply log10 transformation to the field `V`
        V_plot = log10.(V)
        # Color bounds based on the log10-transformed data
        cmin, cmax = minimum(V_plot), maximum(V_plot)
        # Define ticks and labels for the colorbar
        ticks = range(cmin, cmax, length=num_ticks)
        tick_labels = [10^t for t in ticks]  # Exponentially spaced tick labels
        tick_format = "%.1e"
    else
        # Use the original linear values for the field `V`
        V_plot = V
        # Color bounds based on the linear data
        cmin, cmax = minimum(V), maximum(V)
        # Define ticks for the colorbar
        ticks = range(cmin, cmax, length=num_ticks)
        tick_labels = ticks
        tick_format = "%.1e"
    end

    # Create a colorbar for the field `V`
    colormap = :heat
    colorbar = Colorbar(fig[1, 2], limits = (cmin, cmax), colormap = colormap, ticks=ticks)

    # Create a custom layout for sliders with arrow buttons
    slider_layout = GridLayout(fig[2, 1])
    
    # Define slider parameters
    slider_params = [
        (label = "x axis", range = 1:length(x)),
        (label = "y axis", range = 1:length(y)),
        (label = "z axis", range = 1:length(z))
    ]
    
    sliders = []
    left_buttons = []
    right_buttons = []
    value_labels = []
    
    for (i, (label, range)) in enumerate(slider_params)
        # Create left arrow button
        left_btn = Button(slider_layout[i, 1], label = "◀", width = 30)
        push!(left_buttons, left_btn)
        
        # Create slider
        slider = Slider(slider_layout[i, 2:4], range = range, startvalue = range[1])
        push!(sliders, slider)
        
        # Create right arrow button  
        right_btn = Button(slider_layout[i, 5], label = "▶", width = 30)
        push!(right_buttons, right_btn)
        
        # Add label
        Label(slider_layout[i, 6], label, tellwidth = false, halign = :left)
        
        # Add value label that shows current slider value
        value_label = Label(slider_layout[i, 7], string(slider.value[]), tellwidth = false, halign = :center, width = 50)
        push!(value_labels, value_label)
        
        # Update value label when slider changes
        on(slider.value) do val
            value_label.text = string(val)
        end
        
        # Connect arrow buttons to slider
        on(left_btn.clicks) do n
            current_val = slider.value[]
            current_idx = findfirst(x -> x == current_val, range)
            if current_idx !== nothing && current_idx > 1
                set_close_to!(slider, range[current_idx - 1])
            end
        end
        
        on(right_btn.clicks) do n
            current_val = slider.value[]
            current_idx = findfirst(x -> x == current_val, range)
            if current_idx !== nothing && current_idx < length(range)
                set_close_to!(slider, range[current_idx + 1])
            end
        end
    end
    
    # Extract individual sliders for easier reference
    sl_yz, sl_xz, sl_xy = sliders

    # Plot the slices with either log or linear data
    plt = volumeslices!(ax, x, y, z, V_plot, colormap=colormap, colorrange=(cmin, cmax))

    # Set the aspect ratio to be cubic
    # set_aspect_ratio!(ax, 1.0) 

    # Connect sliders to `volumeslices` update methods
    on(sl_yz.value) do v; plt[:update_yz][](v) end
    on(sl_xz.value) do v; plt[:update_xz][](v) end
    on(sl_xy.value) do v; plt[:update_xy][](v) end

    set_close_to!(sl_yz, .5length(x))
    set_close_to!(sl_xz, .5length(y))
    set_close_to!(sl_xy, .5length(z))

    # Add toggles to show/hide heatmaps
    hmaps = [plt[Symbol(:heatmap_, s)][] for s ∈ (:yz, :xz, :xy)]
    toggles = [Toggle(slider_layout[i, 8], active = true) for i ∈ 1:length(hmaps)]

    # Use on() to listen to toggle changes and set visibility directly
    for (hmap, toggle) in zip(hmaps, toggles)
        on(toggle.active) do is_active
            hmap.visible = is_active
        end
    end

    fig
end

function plot_graphs(material_indices, grid_params, Qel, E, V, filename=nothing)
    
    # Plot the voxel grid
    fig1 = plot_slices(material_indices, grid_params, num_ticks = maximum(material_indices), title = "Media Distribution")
    fig2 = plot_slices(Qel, grid_params, title = "Distribution of Qel [W/m^3]")
    fig3 = plot_slices(E, grid_params, title = "Distribution of E [V/m]")
    fig4 = plot_slices(V, Config.GridParams(
        grid_params.lx, 
        grid_params.ly, 
        grid_params.lz, 
        grid_params.nx + 1, 
        grid_params.ny + 1, 
        grid_params.nz + 1), 
        title = "Distribution of V [V]"
    )

    window1 = display(GLMakie.Screen(), fig1)
    window2 = display(GLMakie.Screen(), fig2)
    window3 = display(GLMakie.Screen(), fig3)
    window4 = display(GLMakie.Screen(), fig4)

    # Save figures if filename is provided
    if filename !== nothing
        save("$(filename)_a.png", fig1)
        save("$(filename)_b.png", fig2)
        save("$(filename)_c.png", fig3)
        save("$(filename)_d.png", fig4)
    end

end

end