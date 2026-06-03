module RFSolver

using Ferrite, SparseArrays
using LinearAlgebra
using Statistics
using IterativeSolvers
using ..Config
using ..GridSetup

export solve_rf, calculate_fields, convert_V

function cell_index_to_xyz(cell_number, nx_max, ny_max)
    cell_number -= 1
    x = (cell_number % (nx_max * ny_max)) % nx_max + 1
    y = div((cell_number % (nx_max * ny_max)), nx_max) + 1
    z = div(cell_number, nx_max * ny_max) + 1
    return x, y, z
end

# Function to solve the RF problem in 3D using the finite element method
function solve_rf(grid, rf_params::Config.RFParams, grid_params::Config.GridParams, boundary_conditions)

    σ = rf_params.sigma      # Electrical conductivity matrix [S/m]
    ε_im = rf_params.epsilon_im    # Permittivity matrix [F/m]
    frequency = rf_params.ω          # Frequency of the RF signal [Hz]
    
    # Define the basis and quadrature for the FEM. 
    basis = Lagrange{Ferrite.RefHexahedron, 1}()
    quad = QuadratureRule{Ferrite.RefHexahedron}(2)
    cellvalues = CellValues(quad, basis);
    
    # Initialize the degree of freedom (DoF) handler
    dh = DofHandler(grid)
    add!(dh, :V, basis)  # Add a scalar field 'V' to the DoF handler (for potential)
    close!(dh)

    # Create a sparse stiffness matrix 'K'
    K = allocate_matrix(dh)
    
    # Set up the constraint handler for Dirichlet boundary conditions
    ch = ConstraintHandler(dh)

    # Apply boundary conditions
    for (boundary_name, condition) in boundary_conditions
        add!(ch, Dirichlet(:V, getfacetset(grid, boundary_name), condition))
    end

    close!(ch)

    # Function to assemble the element stiffness matrix and force vector
    # This function computes contributions from the local element to the global matrix
    function assemble_element!(Ke, fe, cellvalues, σ, ε, frequency, cell_number, grid_params)
        n_basefuncs = getnbasefunctions(cellvalues)  # Number of basis functions per element
        fill!(Ke, 0)  # Zero out the element stiffness matrix
        fill!(fe, 0)  # Zero out the element force vector
        
        # Extract the 3D coordinates of the current element

        # Calculate (nx, ny, nz)
        nx_max = grid_params.nx
        ny_max = grid_params.ny

        x, y, z = cell_index_to_xyz(cell_number, nx_max, ny_max)
        
        # Ensure x, y, z are valid indices
        if x > size(σ, 1) || y > size(σ, 2) || z > size(σ, 3)
            throw(ArgumentError("Cell coordinates out of bounds for conductivity or permittivity matrix"))
        end

        conductivity = σ[x, y, z] + frequency * ε_im[x, y, z]
        
        # Loop over all quadrature points for numerical integration
        for q_point in 1:getnquadpoints(cellvalues)
            dΩ = getdetJdV(cellvalues, q_point)  # Jacobian determinant for volume element
            for i in 1:n_basefuncs
                δV  = shape_value(cellvalues, q_point, i)    # Basis function value at q-point
                ∇δV = shape_gradient(cellvalues, q_point, i) # Gradient of the basis function
                # Add contribution to fe
                fe[i] += 0 # δV * dΩ  # Element force vector contribution
                
                # Loop over the basis functions to compute the stiffness matrix contributions
                for j in 1:n_basefuncs
                    ∇V = shape_gradient(cellvalues, q_point, j)
                    # Use correct 3D indexing to access σ and ε
                    # σ[x, y, z] + frequency * ε_im[x, y, z]
                    Ke[i, j] += conductivity * (∇δV ⋅ ∇V) * dΩ
                end
            end
        end
        return Ke, fe
    end
    
    # Function to assemble the global stiffness matrix and force vector
    # This assembles contributions from all elements in the grid
    function assemble_global(cellvalues, K, dh, σ, ε_im, frequency, grid_params)
        # Allocate the element stiffness matrix and element force vector
        n_basefuncs = getnbasefunctions(cellvalues)
        Ke = zeros(n_basefuncs, n_basefuncs)  # Element stiffness matrix
        fe = zeros(n_basefuncs)               # Element force vector
        
        # Allocate global force vector 'f'
        f = zeros(ndofs(dh))  # Global force vector
        # Create an assembler to handle assembly into global matrix/vector
        assembler = start_assemble(K, f)
        
        
        cell_number = 0
        # Loop over all cells (elements) in the grid
        for cell in CellIterator(dh)
            cell_number += 1
            # Reinitialize cell values for this cell
            reinit!(cellvalues, cell)
            # Compute element contributions (Ke and fe)
            assemble_element!(Ke, fe, cellvalues, σ, ε_im, frequency, cell_number, grid_params)
            # Assemble the element contributions into the global matrix and vector
            assemble!(assembler, celldofs(cell), Ke, fe)
        end
        return K, f
    end

    # Assemble the global matrix and solve the system
    K, f = assemble_global(cellvalues, K, dh, σ, ε_im, frequency, grid_params)

    # Apply the Dirichlet boundary conditions
    apply!(K, f, ch)
    
    # Save the vector f to a file
    # open("output_f.txt", "w") do file
    #     for value in f
    #         write(file, "$value\n")
    #     end
    # end

    # Solve the linear system K * V = f for the potential 'V'
    # V = K \ f;
    V = cg(K, f)
    # V = minres(K, f, reltol=1e-12)
    # V = gmres(K, f)
    # Explicitly make sure bcs are correct
    apply!(V, ch)

    return V, dh, cellvalues
end

# Function to calculate the electrical field (E) and power dissipation (Qel) from the potential

function calculate_fields(cellvalues::CellValues, dh::DofHandler, V::AbstractVector{T}, sigma, grid_params) where T
    n = getnbasefunctions(cellvalues)
    cell_dofs = zeros(Int, n)
    nqp = getnquadpoints(cellvalues)

    Q_el = similar(sigma)
    E_new = similar(sigma)

    for (cell_num, cell) in enumerate(CellIterator(dh))
        celldofs!(cell_dofs, dh, cell_num)
        aᵉ = V[cell_dofs]
        reinit!(cellvalues, cell)

        x, y, z = cell_index_to_xyz(cell_num, grid_params.nx, grid_params.ny)

        value = zero(T)
        for q_point in 1:nqp
            g = function_gradient(cellvalues, q_point, aᵉ)
            value += g[1]^2 + g[2]^2 + g[3]^2
        end
        value /= nqp

        Q_el[x, y, z] = (1/2) * sigma[x, y, z] * value
        E_new[x, y, z] = sqrt(value)
    end

    return Q_el, E_new
end


function convert_V(V, grid_params)
    nx_max = grid_params.nx + 1
    ny_max = grid_params.ny + 1
    nz_max = grid_params.nz + 1
    return reshape(V, nx_max, ny_max, nz_max)
end

end # module RFSolver_CPU

