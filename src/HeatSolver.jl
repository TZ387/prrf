# HeatSolver.jl
#
# Solves the classic (non-perfused) heat equation on the structured voxel grid:
#
#   VHC(x) · ∂T/∂t = ∇·(k(x) ∇T) + Q(x)
#
# Spatial discretisation: cell-centred finite differences on the uniform
# Cartesian grid that GridSetup produces.
#
# Time integration: explicit forward Euler, with Δt chosen automatically to
# satisfy the von-Neumann stability criterion:
#
#   Δt ≤ min_over_cells( VHC / (2·k · (1/Δx² + 1/Δy² + 1/Δz²)) )
#
# A safety factor of 0.45 (< 0.5) is applied.
#
# Boundary conditions: zero-flux (Neumann) on all six faces.
#
# ── CPU parallelisation ───────────────────────────────────────────────────────
#
# The solver is called from inside a Threads.@spawn task in TimelapseCreation.
# Threads.@threads is NOT safe when nested inside @spawn — it deadlocks.
#
# Instead, parallelism is achieved by manually spawning one Task per z-slice
# chunk using Threads.@spawn inside heat_laplacian! and euler_update!, then
# waiting on all of them.  This is safe to call from any task or thread.
#
# Worker count is capped at physical CPU cores (Sys.CPU_THREADS ÷ 2).
# Hyperthreads share the FPU and give no benefit for dense FP arithmetic.
# The cap is also bounded by Threads.nthreads() — the Julia threadpool size.
#
# ── GPU pathway (future) ──────────────────────────────────────────────────────
#
# The two hot kernels (heat_laplacian!, euler_update!) are dispatched through
# an AbstractHeatBackend type.  Adding GPU support later means:
#   1. Define  struct CUDABackend <: AbstractHeatBackend end
#   2. Add CuArray methods for heat_laplacian! and euler_update!
#   3. Pass  backend = CUDABackend()  to solve_heat_phase
# Nothing else in this file or in TimelapseCreation needs to change.

module HeatSolver

using ..Config

export compute_stable_dt, solve_heat_phase, cpu_backend

# ── Backend abstraction ───────────────────────────────────────────────────────

abstract type AbstractHeatBackend end

struct CPUBackend <: AbstractHeatBackend
    n_workers :: Int
end

"""
    cpu_backend() -> CPUBackend

Choose worker count at runtime: min(julia_threads, physical_cores).
Physical cores = Sys.CPU_THREADS ÷ 2 (excludes hyperthreads, which share
the FPU and provide no benefit for dense floating-point work).
"""
function cpu_backend()
    physical = max(1, Sys.CPU_THREADS ÷ 2)
    workers  = min(Threads.nthreads(), physical)
    return CPUBackend(workers)
end

# ── Stability ─────────────────────────────────────────────────────────────────

"""
    compute_stable_dt(heat_params, grid_params; safety=0.45) -> Float64

Return the largest Δt satisfying the explicit-Euler von-Neumann stability
condition. Serial loop — called once per phase, not per step.
"""
function compute_stable_dt(heat_params::Config.HeatParams,
                           grid_params::Config.GridParams;
                           safety::Float64 = 0.45)
    dx = grid_params.lx / grid_params.nx
    dy = grid_params.ly / grid_params.ny
    dz = grid_params.lz / grid_params.nz

    inv_sum = 1.0/dx^2 + 1.0/dy^2 + 1.0/dz^2

    k   = heat_params.k
    VHC = heat_params.VHC

    dt_min = Inf
    for i in eachindex(k)
        ki = k[i]
        ki == 0.0 && continue
        dt_cell = VHC[i] / (2.0 * ki * inv_sum)
        dt_min  = min(dt_min, dt_cell)
    end

    return safety * dt_min
end

# ── Helpers ───────────────────────────────────────────────────────────────────

@inline function harmonic_mean(a::Float64, b::Float64)
    s = a + b
    s == 0.0 && return 0.0
    return 2.0 * a * b / s
end

# Compute the laplacian kernel for a single z-slice l.
# Extracted so it can be called from any task without nesting @threads.
@inline function _laplacian_slice!(dT, T, k, nx, ny, l,
                                   idx2::Float64, idy2::Float64, idz2::Float64)
    nz = size(T, 3)
    @inbounds for j in 1:ny, i in 1:nx
        Tc = T[i, j, l]
        kc = k[i, j, l]

        flux_e = i < nx ? harmonic_mean(kc, k[i+1, j, l]) * (T[i+1, j, l] - Tc) : 0.0
        flux_w = i > 1  ? harmonic_mean(k[i-1, j, l], kc) * (T[i-1, j, l] - Tc) : 0.0
        flux_n = j < ny ? harmonic_mean(kc, k[i, j+1, l]) * (T[i, j+1, l] - Tc) : 0.0
        flux_s = j > 1  ? harmonic_mean(k[i, j-1, l], kc) * (T[i, j-1, l] - Tc) : 0.0
        flux_t = l < nz ? harmonic_mean(kc, k[i, j, l+1]) * (T[i, j, l+1] - Tc) : 0.0
        flux_b = l > 1  ? harmonic_mean(k[i, j, l-1], kc) * (T[i, j, l-1] - Tc) : 0.0

        dT[i, j, l] = (flux_e + flux_w) * idx2 +
                      (flux_n + flux_s) * idy2 +
                      (flux_t + flux_b) * idz2
    end
end

# ── Parallel kernels (CPU) ────────────────────────────────────────────────────
#
# We partition the nz z-slices into n_workers contiguous chunks and spawn one
# Task per chunk.  Threads.@spawn is safe to call from inside another spawned
# task (unlike @threads).  wait.(tasks) blocks until all chunks finish.
# No locks needed: each task writes to a disjoint range of dT[:,:,l_start:l_end].

function heat_laplacian!(dT::Array{Float64,3},
                          T::Array{Float64,3},
                          k::Array{Float64,3},
                          dx::Float64, dy::Float64, dz::Float64,
                          backend::CPUBackend)
    nx, ny, nz = size(T)
    idx2 = 1.0 / dx^2
    idy2 = 1.0 / dy^2
    idz2 = 1.0 / dz^2
    nw   = backend.n_workers

    # Partition z-slices into nw chunks
    chunk = max(1, nz ÷ nw)
    tasks = Vector{Task}(undef, nw)

    for w in 1:nw
        l_start = (w - 1) * chunk + 1
        l_end   = (w == nw) ? nz : w * chunk   # last worker takes remainder
        tasks[w] = Threads.@spawn begin
            for l in l_start:l_end
                _laplacian_slice!(dT, T, k, nx, ny, l, idx2, idy2, idz2)
            end
        end
    end

    wait.(tasks)
end

function euler_update!(T::Array{Float64,3},
                       dT::Array{Float64,3},
                       Q::Array{Float64,3},
                       VHC::Array{Float64,3},
                       Δt::Float64,
                       backend::CPUBackend)
    n  = length(T)
    nw = backend.n_workers

    chunk = max(1, n ÷ nw)
    tasks = Vector{Task}(undef, nw)

    for w in 1:nw
        i_start = (w - 1) * chunk + 1
        i_end   = (w == nw) ? n : w * chunk
        tasks[w] = Threads.@spawn begin
            @inbounds for idx in i_start:i_end
                T[idx] += Δt / VHC[idx] * (dT[idx] + Q[idx])
            end
        end
    end

    wait.(tasks)
end

# ── Phase solver ──────────────────────────────────────────────────────────────

"""
    solve_heat_phase(T_in, Q_source, heat_params, grid_params, t_phase;
                     n_update=0, update_cb=nothing, backend=cpu_backend()) -> T

Advance the temperature field over `t_phase` seconds with source `Q_source`.

`backend` selects the compute backend.  Defaults to `cpu_backend()`, which
reads Sys.CPU_THREADS at call time and logs the worker count once per phase.
Pass a future GPU backend here — nothing else in this function needs to change.
"""
function solve_heat_phase(T_in::Array{Float64,3},
                          Q_source::Array{Float64,3},
                          heat_params::Config.HeatParams,
                          grid_params::Config.GridParams,
                          t_phase::Float64;
                          n_update::Int  = 0,
                          update_cb      = nothing,
                          backend        = cpu_backend())

    T  = copy(T_in)
    dT = similar(T)

    dx = grid_params.lx / grid_params.nx
    dy = grid_params.ly / grid_params.ny
    dz = grid_params.lz / grid_params.nz

    k   = heat_params.k
    VHC = heat_params.VHC

    Δt        = compute_stable_dt(heat_params, grid_params)
    num_steps = ceil(Int, t_phase / Δt)
    Δt        = t_phase / num_steps

    @info "Heat phase: t_phase=$(t_phase) s, " *
          "Δt=$(round(Δt; sigdigits=4)) s, " *
          "steps=$num_steps, " *
          "workers=$(backend.n_workers)"

    update_steps = Set{Int}()
    if n_update > 0 && update_cb !== nothing
        for i in 1:n_update
            push!(update_steps, round(Int, i * num_steps / n_update))
        end
    end

    t = 0.0
    for step in 1:num_steps
        heat_laplacian!(dT, T, k, dx, dy, dz, backend)
        euler_update!(T, dT, Q_source, VHC, Δt, backend)
        t += Δt
        if step in update_steps
            update_cb(T, t)
        end
    end

    return T
end

end # module HeatSolver