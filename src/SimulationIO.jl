module SimulationIO

using HDF5
using ..Config

export save_simulation, load_simulation

function save_simulation(filepath::String, grid_params::Config.GridParams, material_indices, Qel, E_new, V_new)
    h5open(filepath, "w") do f
        # Grid parameters as attributes on a dedicated group
        g = create_group(f, "grid_params")
        g["lx"] = grid_params.lx
        g["ly"] = grid_params.ly
        g["lz"] = grid_params.lz
        g["nx"] = grid_params.nx
        g["ny"] = grid_params.ny
        g["nz"] = grid_params.nz

        f["material_indices"] = material_indices
        f["Qel"]              = Qel
        f["E_new"]            = E_new
        f["V_new"]            = V_new
    end
end

function load_simulation(filepath::String)
    h5open(filepath, "r") do f
        g = f["grid_params"]
        grid_params = Config.GridParams(
            read(g["lx"]), read(g["ly"]), read(g["lz"]),
            read(g["nx"]), read(g["ny"]), read(g["nz"])
        )

        material_indices = read(f["material_indices"])
        Qel              = read(f["Qel"])
        E_new            = read(f["E_new"])
        V_new            = read(f["V_new"])

        return grid_params, material_indices, Qel, E_new, V_new
    end
end

end # module SimulationIO
