# TimelapseCreation.jl
#
# Live-updating heat simulation visualisation using GLMakie observables.
#
# The simulation follows the arbitrary schedule defined in HeatParams.schedule —
# an ordered list of (:on/:off, duration) phases, e.g.:
#
#   [(:on, 30.0), (:off, 60.0), (:on, 15.0)]
#
# The cross-section plot (XZ mid-plane slice of T) is refreshed exactly
# n_update times per phase.  Plot updates happen inside the time-stepping loop
# via a callback — no frames are captured; the window just updates live.

module TimelapseCreation

using GLMakie
using Printf
using ..Config
using ..BioheatSolver

export run_heat_timelapse

# ── Helper: mid-plane slice index ────────────────────────────────────────────

mid(n) = max(1, n ÷ 2)

# ── Main entry point ──────────────────────────────────────────────────────────

"""
    run_heat_timelapse(Qel, heat_params, grid_params)

Run the heat simulation according to `heat_params.schedule` and display a
live-updating XZ cross-section plot of the temperature field.

The RF source `Qel` is an (nx, ny, nz) array of volumetric power density
[W/m³] computed by the RF solver.  It is used during `:on` phases; `:off`
phases use a zero source.

The plot is refreshed `n_update` times per phase.
"""
function run_heat_timelapse(Qel::Array{Float64,3},
                            heat_params::Config.HeatParams,
                            grid_params::Config.GridParams)

    nx, ny, nz = grid_params.nx, grid_params.ny, grid_params.nz
    jmid = mid(ny)
    Qzero = zeros(Float64, nx, ny, nz)

    # ── Initial temperature field ─────────────────────────────────────────────
    T = fill(heat_params.T_initial, nx, ny, nz)

    # ── Makie scene setup ─────────────────────────────────────────────────────
    slice_obs  = Observable(T[:, jmid, :])
    title_obs  = Observable("t = 0.00 s")
    crange_obs = Observable((heat_params.T_initial, heat_params.T_initial + 1.0))

    fig = Figure(size = (800, 600))
    ax  = Axis(fig[1, 1];
               title  = title_obs,
               xlabel = "x index",
               ylabel = "z index")

    hm = heatmap!(ax, slice_obs;
                  colormap   = :thermal,
                  colorrange = crange_obs)
    Colorbar(fig[1, 2], hm; label = "Temperature [°C]")

    display(fig)

    # ── Callback factory ──────────────────────────────────────────────────────
    # Returns a closure over (phase_label, t_offset) that updates the
    # observables and triggers a repaint when called with (T, t_local).

    function make_callback(phase_label::String, t_offset::Float64)
        return function(T_current::Array{Float64,3}, t_local::Float64)
            Tslice = T_current[:, jmid, :]
            Tmin   = minimum(T_current)
            Tmax   = maximum(T_current)
            if Tmax ≈ Tmin
                Tmax = Tmin + 1.0
            end
            slice_obs[]  = Tslice
            crange_obs[] = (Tmin, Tmax)
            title_obs[]  = @sprintf("t = %.2f s  [%s]", t_offset + t_local, phase_label)
            sleep(0.0)  # yield to the event loop so GLMakie repaints
        end
    end

    # ── Schedule loop ─────────────────────────────────────────────────────────
    t_elapsed = 0.0   # wall-clock simulation time across all phases

    for (phase_idx, (state, duration)) in enumerate(heat_params.schedule)

        state in (:on, :off) || error(
            "Unknown phase state $(repr(state)) in schedule entry $phase_idx. " *
            "Expected :on or :off.")

        label  = state == :on ? "heating" : "cooling"
        Q_src  = state == :on ? Qel : Qzero

        @info "Phase $phase_idx/$( length(heat_params.schedule)): " *
              "$label for $(duration) s  ($(heat_params.n_update) plot updates)"

        cb = make_callback(label, t_elapsed)
        T  = BioheatSolver.solve_heat_phase(
                 T, Q_src, heat_params, grid_params, duration;
                 n_update  = heat_params.n_update,
                 update_cb = cb)

        @info "Phase $phase_idx complete.  Peak T = $(round(maximum(T); digits=2)) °C"

        t_elapsed += duration
    end

    return T
end

end # module TimelapseCreation