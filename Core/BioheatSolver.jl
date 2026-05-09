# BioheatSolver.jl
module BioheatSolver

using LinearAlgebra
using Ferrite
using Tensors

export solve_bioheat

function solve_bioheat(grid, T_old, Qel, bioheat_params, grid_params)
    # Define basis and quadrature for FEM
    basis = Lagrange{3, RefCube, 1}()
    quad = QuadratureRule{3, RefCube}(2)
    dh = DofHandler(grid)
    push!(dh, :T, basis)
    close!(dh)

    # Assemble the thermal conductivity matrix K_t
    K_t_sp = create_sparsity_pattern(dh)
    K_t = copy(K_t_sp)
    assembler_t = start_assemble(K_t)

    for elem in 1:length(grid)
        cell = get_cell(grid, elem)
        dofs = get_element_dofs(dh, cell)
        Ke_t = zeros(Float64, ndofs(dh), ndofs(dh))
        
        for qp in 1:num_quadrature_points(quad)
            q_point = quadrature_point(quad, qp)
            detJ = jacobian_determinant(cell, q_point)
            B = gradient_basis_functions(basis, cell, q_point)
            Ke_t += B' * bioheat_params.k * B * detJ * weight(quad, qp)
        end
        
        assemble!(assembler_t, dofs, Ke_t)
    end
    
    K_t = end_assemble(assembler_t)

    # Assemble the mass matrix M
    M_sp = create_sparsity_pattern(dh)
    M = copy(M_sp)
    assembler_m = start_assemble(M)

    for elem in 1:length(grid)
        cell = get_cell(grid, elem)
        dofs = get_element_dofs(dh, cell)
        Me = zeros(Float64, ndofs(dh), ndofs(dh))
        
        for qp in 1(num_quadrature_points(quad))
            q_point = quadrature_point(quad, qp)
            detJ = jacobian_determinant(cell, q_point)
            N = basis_functions(basis, cell, q_point)
            Me += N' * bioheat_params.rho * bioheat_params.c * N * detJ * weight(quad, qp)
        end
        
        assemble!(assembler_m, dofs, Me)
    end
    
    M = end_assemble(assembler_m)

    # Time-stepping loop
    T = zeros(Float64, ndofs(dh))
    for step in 1:bioheat_params.num_steps
        # Assemble the source term vector F
        F = zeros(Float64, ndofs(dh))
        for elem in 1:length(grid)
            cell = get_cell(grid, elem)
            dofs = get_element_dofs(dh, cell)
            for qp in 1(num_quadrature_points(quad))
                q_point = quadrature_point(quad, qp)
                N = basis_functions(basis, cell, q_point)
                F[dofs] += (bioheat_params.Qmet + Qel[elem] + bioheat_params.rho_b * bioheat_params.cb * bioheat_params.ωb * (bioheat_params.Ta - T_old[dofs])) * N * weight(quad, qp)
            end
        end
        
        # Solve the system
        A = M / bioheat_params.Δt + K_t
        b = M * T_old / bioheat_params.Δt + F
        T = A \ b
        
        # Update temperature for the next time step
        T_old = T
    end

    return T
end

end # module BioheatSolver
