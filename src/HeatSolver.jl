# HeatSolver.jl
#
# Solves the classic (non-perfused) heat equation on the structured voxel grid:
#
#   VHC(x) · ∂T/∂t = ∇·(k(x) ∇T) + Q(x)
#
# Spatial discretisation: cell-centred finite differences on the uniform
# Cartesian grid that GridSetup produces.  This is a natural fit for the
# existing (nx, ny, nz) cell arrays and avoids the complexity of a full FEM
# assembly.
#
# Time integration: explicit forward Euler, with Δt chosen automatically to
# satisfy the von-Neumann stability criterion:
#
#   Δt ≤ min_over_cells( VHC / (2·k · (1/Δx² + 1/Δy² + 1/Δz²)) )
#
# A safety factor of 0.45 (< 0.5) is applied.
#
# Boundary conditions: zero-flux (Neumann) on all six faces, implemented via
# ghost-cell reflection (one-sided differences at the boundary).

module HeatSolver

using ..Config

export compute_stable_dt, solve_heat_phase

# ── Stability ────────────────────────────────────────────────────────────────

"""
    compute_stable_dt(heat_params, grid_params; safety=0.45) -> Float64

Return the largest Δt that satisfies the explicit-Euler von-Neumann stability
condition for the heterogeneous heat equation on a uniform Cartesian grid.

The criterion per cell is  Δt ≤ VHC / (2k · Σ 1/Δξ²),  and we take the
global minimum over all cells.  `safety < 0.5` provides a margin.
"""
function compute_stable_dt(heat_params::Config.HeatParams,
                           grid_params::Config.GridParams;
                           safety::Float64 = 0.45)
    dx = grid_params.lx / grid_params.nx
    dy = grid_params.ly / grid_params.ny
    dz = grid_params.lz / grid_params.nz

    inv_sum = 1.0/dx^2 + 1.0/dy^2 + 1.0/dz^2   # same for every cell (uniform grid)

    k   = heat_params.k
    VHC = heat_params.VHC

    # Avoid division by zero for cells with k = 0 (should not happen physically,
    # but guard anyway).
    dt_min = Inf
    for i in eachindex(k)
        ki = k[i]
        ki == 0.0 && continue
        dt_cell = VHC[i] / (2.0 * ki * inv_sum)
        dt_min  = min(dt_min, dt_cell)
    end

    return safety * dt_min
end

# ── Laplacian helper ─────────────────────────────────────────────────────────

# Central-difference divergence of the heat flux  ∇·(k ∇T)  evaluated at cell
# (i,j,l) using the harmonic mean conductivity at the cell interface:
#
#   k_{i+½} = 2 k_i k_{i+1} / (k_i + k_{i+1})
#
# Zero-flux (Neumann) BCs are imposed via one-sided differences:
# the ghost value beyond the boundary equals the boundary cell value, so the
# interface conductance effectively becomes zero there (∂T/∂n = 0).

@inline function harmonic_mean(a::Float64, b::Float64)
    s = a + b
    s == 0.0 && return 0.0
    return 2.0 * a * b / s
end

function heat_laplacian!(dT::Array{Float64,3},
                         T::Array{Float64,3},
                         k::Array{Float64,3},
                         dx::Float64, dy::Float64, dz::Float64)
    nx, ny, nz = size(T)

    @inbounds for l in 1:nz, j in 1:ny, i in 1:nx
        Tc = T[i, j, l]
        kc = k[i, j, l]

        # ── x-direction ──────────────────────────────────────────────────────
        # east interface
        if i < nx
            ke  = harmonic_mean(kc, k[i+1, j, l])
            flux_e = ke * (T[i+1, j, l] - Tc)
        else
            flux_e = 0.0  # zero-flux BC
        end
        # west interface
        if i > 1
            kw  = harmonic_mean(k[i-1, j, l], kc)
            flux_w = kw * (T[i-1, j, l] - Tc)
        else
            flux_w = 0.0
        end

        # ── y-direction ──────────────────────────────────────────────────────
        if j < ny
            kn  = harmonic_mean(kc, k[i, j+1, l])
            flux_n = kn * (T[i, j+1, l] - Tc)
        else
            flux_n = 0.0
        end
        if j > 1
            ks  = harmonic_mean(k[i, j-1, l], kc)
            flux_s = ks * (T[i, j-1, l] - Tc)
        else
            flux_s = 0.0
        end

        # ── z-direction ──────────────────────────────────────────────────────
        if l < nz
            kt  = harmonic_mean(kc, k[i, j, l+1])
            flux_t = kt * (T[i, j, l+1] - Tc)
        else
            flux_t = 0.0
        end
        if l > 1
            kb  = harmonic_mean(k[i, j, l-1], kc)
            flux_b = kb * (T[i, j, l-1] - Tc)
        else
            flux_b = 0.0
        end

        dT[i, j, l] = (flux_e + flux_w) / dx^2 +
                      (flux_n + flux_s) / dy^2 +
                      (flux_t + flux_b) / dz^2
    end
end

# ── Phase solver ─────────────────────────────────────────────────────────────

"""
    solve_heat_phase(T, Qel_source, heat_params, grid_params, t_phase;
                     n_update=0, update_cb=nothing) -> T

Advance the temperature field `T` (mutated in-place, a copy is made internally)
over `t_phase` seconds using the heat equation with source `Qel_source`.

`n_update` controls how many times the optional callback `update_cb(T, t)` is
called during the phase.  Set `n_update = 0` to suppress all callbacks.

Returns the updated temperature array.
"""
function solve_heat_phase(T_in::Array{Float64,3},
                          Qel_source::Array{Float64,3},
                          heat_params::Config.HeatParams,
                          grid_params::Config.GridParams,
                          t_phase::Float64;
                          n_update::Int = 0,
                          update_cb = nothing)

    T   = copy(T_in)
    dT  = similar(T)

    dx  = grid_params.lx / grid_params.nx
    dy  = grid_params.ly / grid_params.ny
    dz  = grid_params.lz / grid_params.nz

    k   = heat_params.k
    VHC = heat_params.VHC

    Δt  = compute_stable_dt(heat_params, grid_params)

    # Total number of time steps for this phase
    num_steps = ceil(Int, t_phase / Δt)
    # Recompute exact Δt so the phase duration is hit exactly
    Δt = t_phase / num_steps

    @info "Heat phase: t_phase=$(t_phase) s, Δt=$(round(Δt; sigdigits=4)) s, steps=$num_steps"

    # Decide at which step indices we call the update callback.
    # We want n_update evenly-spaced snapshots across the phase.
    update_steps = Set{Int}()
    if n_update > 0 && update_cb !== nothing
        for i in 1:n_update
            push!(update_steps, round(Int, i * num_steps / n_update))
        end
    end

    t = 0.0
    for step in 1:num_steps

        # Compute ∇·(k ∇T) into dT
        heat_laplacian!(dT, T, k, dx, dy, dz)

        # Forward Euler update:  T += Δt/VHC * (∇·(k∇T) + Q)
        @inbounds for idx in eachindex(T)
            T[idx] += Δt / VHC[idx] * (dT[idx] + Qel_source[idx])
        end

        t += Δt

        if step in update_steps
            update_cb(T, t)
        end
    end

    return T
end

end # module HeatSolver