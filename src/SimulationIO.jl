module SimulationIO

using HDF5
using ..Config

export save_simulation, load_simulation

function save_simulation(filepath::String,
                         grid_params::Config.GridParams,
                         material_indices;
                         kwargs...)
    h5open(filepath, "w") do f
        _save_grid_params(f, grid_params)
        f["material_indices"] = material_indices

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

function load_simulation(filepath::String)
    h5open(filepath, "r") do f
        grid_params      = _load_grid_params(f)
        material_indices = read(f["material_indices"])

        data = if haskey(f, "data")
            _load_group(f["data"])
        else
            NamedTuple()
        end

        return grid_params, material_indices, data
    end
end

function _save_grid_params(fid, gp::Config.GridParams)
    g = create_group(fid, "grid_params")
    for fname in (:lx, :ly, :lz, :nx, :ny, :nz)
        g[string(fname)] = getfield(gp, fname)
    end
end

function _load_grid_params(fid)
    g = fid["grid_params"]
    Config.GridParams(
        read(g["lx"]), read(g["ly"]), read(g["lz"]),
        read(g["nx"]), read(g["ny"]), read(g["nz"])
    )
end

function _write_field(group, name::String, v)
    if v isa AbstractArray{<:Complex}
        sg = create_group(group, name)
        attrs(sg)["complex"] = 1
        chunk = _autochunk(size(v))
        _write_compressed(sg, "real", Float64.(real.(v)), chunk)
        _write_compressed(sg, "imag", Float64.(imag.(v)), chunk)
    elseif v isa AbstractArray{<:Real}
        _write_compressed(group, name, Float64.(v), _autochunk(size(v)))
    elseif v isa Number
        group[name] = v
    elseif v isa AbstractString
        group[name] = v
    else
        @warn "SimulationIO: field '$name' has unsupported element type $(eltype(v)); storing as string."
        group[name] = string(v)
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

function _autochunk(sz::Tuple)
    TARGET = 32_768
    if length(sz) == 1
        return (min(sz[1], TARGET),)
    elseif length(sz) == 2
        s = round(Int, sqrt(TARGET))
        return (min(sz[1], s), min(sz[2], s))
    elseif length(sz) == 3
        s = round(Int, cbrt(TARGET))
        return (min(sz[1], s), min(sz[2], s), min(sz[3], s))
    else
        return sz
    end
end

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
                push!(pairs_vec, sym => _load_group(obj))
            end
        else
            push!(pairs_vec, sym => read(obj))
        end
    end
    return NamedTuple(pairs_vec)
end

end # module SimulationIO
