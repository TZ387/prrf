# Timelapse.jl
#
# Live-updating heat simulation visualisation using GLMakie observables.
#
# The simulation follows the arbitrary schedule defined in HeatParams.schedule —
# an ordered list of (:on/:off, duration) phases, e.g.:
#
#   [(:on, 30.0), (:off, 60.0), (:on, 15.0)]
#
# Threading model
# ───────────────
# The root cause of VS Code's "not responding" freeze is that the heavy solver
# loop and the GLMakie window event loop both run on the same (main) thread.
# The OS sees the thread busy for seconds at a time and marks the window
# unresponsive.
#
# Fix: the entire heat simulation runs on a worker thread (Threads.@spawn).
# The main thread never blocks on the solver — it sits in a lightweight
# polling loop that (a) drains incoming plot snapshots from a Channel and
# writes them to Observables, and (b) calls sleep(0.05) between polls to
# yield back to the GLMakie event loop.
#
# Why a Channel and not direct Observable writes from the worker?
# GLMakie requires all Observable assignments to happen on the main thread.
# Writing obs[] = ... from a worker thread causes race conditions and
# occasional crashes.  The Channel decouples the two threads cleanly:
# the worker puts a cheap snapshot struct into it; the main thread takes
# snapshots out and does all Observable writes itself.
#
# Layout:
#   worker thread → solve_heat_phase × N → put!(channel, snapshot)
#   main thread   → polling loop { take! + obs[] = ..., sleep(0.05) }

module Timelapse

using GLMakie
using Printf
using Base.Threads
using ..Config
using ..HeatSolver

export run_heat_timelapse

# ── Snapshot carried through the channel ─────────────────────────────────────

struct PlotSnapshot
    slice :: Matrix{Float64}  # copy of T[:, jmid, :], made on worker thread
    t_sim :: Float64          # simulation time [s]
    label :: String
    T_min :: Float64
    T_max :: Float64
end

# ── Helper ────────────────────────────────────────────────────────────────────

mid(n) = max(1, n ÷ 2)

# ── Main entry point ──────────────────────────────────────────────────────────

"""
    run_heat_timelapse(Qel, heat_params, grid_params)

Run the multi-phase heat simulation on a worker thread while displaying a
live-updating XZ cross-section plot on the main thread.

The window stays fully responsive because the main thread never blocks on
the solver — it only services the GLMakie event loop and drains a Channel
of plot snapshots produced by the worker.
"""
function run_heat_timelapse(Qel::Array{Float64,3},
                            heat_params::Config.HeatParams,
                            grid_params::Config.GridParams)

    nx, ny, nz = grid_params.nx, grid_params.ny, grid_params.nz
    jmid  = mid(ny)
    Qzero = zeros(Float64, nx, ny, nz)

    # ── Channel ───────────────────────────────────────────────────────────────
    # Capacity 4: small buffer so the worker is never held up long waiting for
    # the main thread to drain, but also doesn't queue stale frames.
    # Nothing (::Nothing) is the sentinel that signals simulation is done.
    snapshot_ch = Channel{Union{PlotSnapshot, Nothing}}(4)

    # ── Build the Makie scene on the main thread ──────────────────────────────
    T0         = fill(heat_params.T_initial, nx, ny, nz)
    slice_obs  = Observable(T0[:, jmid, :])
    title_obs  = Observable("t = 0.00 s")
    crange_obs = Observable((heat_params.T_initial, heat_params.T_initial + 1.0))

    fig = Figure(size = (800, 600))
    ax  = Axis(fig[1, 1];
               title  = title_obs,
               xlabel = "x index",
               ylabel = "z index")
    hm  = heatmap!(ax, slice_obs;
                   colormap   = :thermal,
                   colorrange = crange_obs)
    Colorbar(fig[1, 2], hm; label = "Temperature [°C]")
    display(fig)

    # ── Worker thread ─────────────────────────────────────────────────────────
    sim_task = Threads.@spawn begin
        T         = copy(T0)
        t_elapsed = 0.0

        for (phase_idx, (state, duration)) in enumerate(heat_params.schedule)

            state in (:on, :off) || error(
                "Unknown phase state $(repr(state)) in schedule entry $phase_idx. " *
                "Expected :on or :off.")

            label = state == :on ? "heating" : "cooling"
            Q_src = state == :on ? Qel : Qzero

            @info "Phase $phase_idx/$(length(heat_params.schedule)): " *
                  "$label for $(duration) s  ($(heat_params.n_update) plot updates)"

            # t_elapsed captured by value for this phase's closure
            t_off = t_elapsed

            cb = (T_curr, t_local) -> begin
                s    = T_curr[:, jmid, :]   # copy slice on worker — safe
                Tmin = minimum(T_curr)
                Tmax = maximum(T_curr)
                Tmax = Tmax ≈ Tmin ? Tmin + 1.0 : Tmax
                put!(snapshot_ch,
                     PlotSnapshot(s, t_off + t_local, label, Tmin, Tmax))
            end

            T = HeatSolver.solve_heat_phase(
                    T, Q_src, heat_params, grid_params, duration;
                    n_update  = heat_params.n_update,
                    update_cb = cb)

            @info "Phase $phase_idx complete.  " *
                  "Peak T = $(round(maximum(T); digits=2)) °C"

            t_elapsed += duration
        end

        put!(snapshot_ch, nothing)  # sentinel: tell main thread we are done
        return T
    end

    # ── Main thread polling loop ──────────────────────────────────────────────
    # Runs until the worker sends the nothing sentinel.
    # sleep(0.05) is the key: it yields to the GLMakie event loop every 50 ms,
    # which is more than enough for a responsive window.
    done = false
    while !done
        # Drain everything currently queued before sleeping
        while isready(snapshot_ch)
            snap = take!(snapshot_ch)
            if snap === nothing
                done = true
                break
            end
            slice_obs[]  = snap.slice
            crange_obs[] = (snap.T_min, snap.T_max)
            title_obs[]  = @sprintf("t = %.2f s  [%s]", snap.t_sim, snap.label)
        end
        sleep(0.05)
    end

    # Re-throw any exception that occurred on the worker, and return T_final.
    return fetch(sim_task)
end

end # module Timelapse