module SimulationIO

using HDF5
using ..Config

export save_simulation, load_simulation

# ============================================================================
# PUBLIC API
# ============================================================================

"""
    save_simulation(filepath, grid_params, material_indices; kwargs...)

Save a simulation to an HDF5 file.

`grid_params` and `material_indices` are always saved in their own fixed
groups.  Any additional data — arrays, scalars, or complex arrays — is passed
as keyword arguments and written into a `data/` group, one dataset per keyword.

# Examples
```julia
# Minimal (original behaviour)
save_simulation("out.h5", grid_params, material_indices;
                Qel=Qel, E_new=E_new, V_new=V_new)

# Extended (Example 6: also save V and E)
save_simulation("out.h5", grid_params, material_indices;
                Qel=Qel, E_new=E_new, V_new=V_new, V=V, E=E)

# Arbitrary fields — add whatever the caller needs
save_simulation("out.h5", grid_params, material_indices;
                temperature=T, pressure=P, my_scalar=42.0)
```

Complex arrays are automatically split into real/imag sub-datasets and
reconstructed transparently on load.  All array datasets are stored with
GZIP compression (level 3) + shuffle filter.
"""
function save_simulation(filepath::String,
                         grid_params::Config.GridParams,
                         material_indices;
                         kwargs...)
    h5open(filepath, "w") do f
        # --- fixed groups ---------------------------------------------------
        _save_grid_params(f, grid_params)
        f["material_indices"] = material_indices

        # --- generic data group --------------------------------------------
        if !isempty(kwargs)
            g = create_group(f, "data")
            attrs(g)["description"] = "Simulation output fields"
            for (k, v) in kwargs
                _write_field(g, string(k), v)
            end
        end
    end
    println("Simulation saved → $filepath")
end

# ============================================================================

"""
    load_simulation(filepath) -> (grid_params, material_indices, data)

Load a simulation saved by `save_simulation`.

Returns:
* `grid_params`       – `Config.GridParams` struct
* `material_indices`  – array
* `data`              – `NamedTuple` containing every field that was passed as
                        a keyword argument to `save_simulation`.  Access fields
                        with `data.Qel`, `data.E_new`, `data.V`, etc.

# Example
```julia
grid_params, material_indices, data = load_simulation("out.h5")
Qel   = data.Qel
E_new = data.E_new
V     = data.V      # only present if it was saved
```
"""
function load_simulation(filepath::String)
    h5open(filepath, "r") do f
        grid_params      = _load_grid_params(f)
        material_indices = read(f["material_indices"])

        data = if haskey(f, "data")
            _load_group(f["data"])
        else
            NamedTuple()   # file has no data group (e.g. old format)
        end

        return grid_params, material_indices, data
    end
end

# ============================================================================
# INTERNAL — GRID PARAMS
# ============================================================================

function _save_grid_params(fid, gp::Config.GridParams)
    g = create_group(fid, "grid_params")
    for fname in (:lx, :ly, :lz, :nx, :ny, :nz)
        g[string(fname)] = getfield(gp, fname)
    end
    attrs(g)["description"] = "GridParams: domain extents and element counts"
end

function _load_grid_params(fid)
    g = fid["grid_params"]
    Config.GridParams(
        read(g["lx"]), read(g["ly"]), read(g["lz"]),
        read(g["nx"]), read(g["ny"]), read(g["nz"])
    )
end

# ============================================================================
# INTERNAL — GENERIC FIELD WRITER
# ============================================================================

"""Write a single named field (scalar, real array, or complex array)."""
function _write_field(group, name::String, v)
    if v isa AbstractArray{<:Complex}
        sg = create_group(group, name)
        attrs(sg)["complex"] = 1
        chunk = _autochunk(size(v))
        _write_compressed(sg, "real", Float64.(real.(v)), chunk)
        _write_compressed(sg, "imag", Float64.(imag.(v)), chunk)
    elseif v isa AbstractArray{<:Real}
        _write_compressed(group, name, Float64.(v), _autochunk(size(v)))
    elseif v isa AbstractArray
        # Fallback for integer arrays etc. — write as-is (no compression helper)
        group[name] = v
    elseif v isa Number
        group[name] = v
    elseif v isa AbstractString
        group[name] = v
    else
        group[name] = string(v)   # last resort
    end
end

function _write_compressed(group, name::String, data::Array{Float64}, chunk)
    ds = create_dataset(
        group, name, datatype(Float64), dataspace(data);
        chunk   = chunk,
        shuffle = (),
        deflate = 3
    )
    write_dataset(ds, datatype(Float64), data)
end

"""Choose a sensible chunk size for compression."""
function _autochunk(sz::Tuple)
    if length(sz) == 1
        return (min(sz[1], 1024),)
    elseif length(sz) == 2
        return (min(sz[1], 64), min(sz[2], 64))
    elseif length(sz) == 3
        return (min(sz[1], 32), min(sz[2], 32), min(sz[3], 32))
    else
        return sz
    end
end

# ============================================================================
# INTERNAL — GENERIC GROUP READER
# ============================================================================

"""Read every dataset/sub-group in an HDF5 group and return a NamedTuple."""
function _load_group(group)
    pairs_vec = Pair{Symbol, Any}[]

    for key in keys(group)
        obj = group[key]
        sym = Symbol(key)

        if obj isa HDF5.Group
            if haskey(attrs(obj), "complex") && read(attrs(obj)["complex"]) == 1
                re  = read(obj["real"])
                im_ = read(obj["imag"])
                push!(pairs_vec, sym => complex.(re, im_))
            else
                push!(pairs_vec, sym => _load_group(obj))   # recurse
            end
        else
            push!(pairs_vec, sym => read(obj))
        end
    end

    return NamedTuple(pairs_vec)
end

end # module SimulationIO