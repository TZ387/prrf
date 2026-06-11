# TimelapseCreation.jl
#
# Live-updating heat simulation visualisation using GLMakie observables.
#
# The simulation is split into two phases driven by BioheatSolver:
#
#   1. "on"  phase (duration t_on):  RF source Qel is active.
#   2. "off" phase (duration t_off): RF source is zero (cooling).
#
# The cross-section plot (XZ mid-plane slice of T) is refreshed exactly
# n_update times during each phase.  Plot updates happen inside the
# time-stepping loop via a callback — no frames are captured; the window
# just updates live.

module TimelapseCreation

using GLMakie
using Printf
using ..Config
using ..BioheatSolver

export run_heat_timelapse

# ── Helper: mid-plane slice indices ─────────────────────────────────────────

mid(n) = max(1, n ÷ 2)

# ── Main entry point ─────────────────────────────────────────────────────────

"""
    run_heat_timelapse(Qel, bioheat_params, grid_params)

Run the two-phase heat simulation (on → off) and display a live-updating
cross-section plot of the temperature field.

The RF source `Qel` is an (nx, ny, nz) array of volumetric power density
[W/m³] computed by the RF solver.

Updates the plot `n_update` times during the heating phase and `n_update`
times during the cooling phase (total 2·n_update frame refreshes).
"""
function run_heat_timelapse(Qel::Array{Float64,3},
                            heat_params::Config.HeatParams,
                            grid_params::Config.GridParams)

    nx, ny, nz = grid_params.nx, grid_params.ny, grid_params.nz
    jmid = mid(ny)   # fixed y-index for the XZ slice

    # ── Initial temperature field ────────────────────────────────────────────
    T = fill(heat_params.T_initial, nx, ny, nz)

    # ── Makie scene setup ────────────────────────────────────────────────────
    # Observable wrapping the 2-D slice shown in the heatmap.
    slice_obs  = Observable(T[:, jmid, :])           # (nx, nz) matrix
    title_obs  = Observable("t = 0.00 s  [heating]")

    fig = Figure(size = (800, 600))
    ax  = Axis(fig[1, 1];
               title  = title_obs,
               xlabel = "x index",
               ylabel = "z index")

    # Color range observable — updated as simulation progresses so the scale
    # always spans the current min/max of the full 3-D field.
    crange_obs = Observable((heat_params.T_initial, heat_params.T_initial + 1.0))

    hm = heatmap!(ax, slice_obs;
                  colormap = :thermal,
                  colorrange = crange_obs)
    Colorbar(fig[1, 2], hm; label = "Temperature [°C]")

    display(fig)

    # ── Callback factory ─────────────────────────────────────────────────────
    # Returns a closure that, when called with (T, t), updates the observables
    # and triggers a GLMakie render.

    function make_callback(phase_label::String)
        return function(T_current::Array{Float64,3}, t::Float64)
            Tslice = T_current[:, jmid, :]
            Tmin   = minimum(T_current)
            Tmax   = maximum(T_current)
            # Avoid degenerate color range if T is uniform
            if Tmax ≈ Tmin
                Tmax = Tmin + 1.0
            end
            slice_obs[]  = Tslice
            crange_obs[] = (Tmin, Tmax)
            title_obs[]  = @sprintf("t = %.2f s  [%s]", t, phase_label)
            # Flush the event queue so the window repaints immediately
            sleep(0.0)
        end
    end

    # ── Phase 1: RF on ───────────────────────────────────────────────────────
    @info "Starting heating phase (t_on = $(heat_params.t_on) s, " *
          "$(heat_params.n_update) plot updates)..."

    cb_on = make_callback("heating")
    T = BioheatSolver.solve_heat_phase(
            T, Qel, heat_params, grid_params, heat_params.t_on;
            n_update = heat_params.n_update,
            update_cb = cb_on)

    @info "Heating phase complete.  Peak T = $(round(maximum(T); digits=2)) °C"

    # ── Phase 2: RF off (cooling) ────────────────────────────────────────────
    @info "Starting cooling phase (t_off = $(heat_params.t_off) s, " *
          "$(heat_params.n_update) plot updates)..."

    Qzero = zeros(Float64, nx, ny, nz)
    cb_off = make_callback("cooling")

    # Offset time so the displayed time continues from t_on
    t_offset = heat_params.t_on
    cb_off_offset = (T_curr, t) -> cb_off(T_curr, t + t_offset)

    T = BioheatSolver.solve_heat_phase(
            T, Qzero, heat_params, grid_params, heat_params.t_off;
            n_update = heat_params.n_update,
            update_cb = cb_off_offset)

    @info "Cooling phase complete.  Final peak T = $(round(maximum(T); digits=2)) °C"

    return T
end

end # module TimelapseCreation