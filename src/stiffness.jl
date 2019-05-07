using LinearAlgebra: norm
using SparseArrays

function local_coord_rotation(u::Vector{Float64}, v::Vector{Float64})
    local L = norm(u - v)
    c = (v - u) ./ L
    R = zeros(Float64, 3, 3)
    R[1,1] = c[1]
    R[1,2] = c[2]
    R[2,1] = -c[2]
    R[2,2] = c[1]
    R[3,3] = 1
    return R
end

"""
local_stiffness_matrix: compute elemental FRAME stiffness matrix
in the local coordinates.

Return: 6 x 6 matrix K_le (local elemental)

Note: only support 2d right now.
"""
function local_stiffness_matrix(L::Float64, E::Float64, A::Float64)
    K_le = zeros(Float64, 2, 2)
    K_le[1,1] = 1.0
    K_le[1,2] = -1.0
    K_le[2,1] = -1.0
    K_le[2,2] = 1.0
    K_le *= E * A / L
    return K_le
end

# TODO: can be optimized using symbolic substitution, only X is changing here
# NOTE: A might be set as ones, if we are using ∑ f_e * l_e formula?
function assemble_global_stiffness_matrix(X::Matrix{Float64}, T::Matrix{Int64},
                                          A::Vector{Float64}, E::Float64,
                                          node_dof::Int, full_node_dof::Int)
    n_elements = size(T,1)
    n_nodes = size(X,1)
    # dimension get from X

    # element -> dof id map
    sys_dof::Int = node_dof * n_nodes
    id_map = zeros(Int, n_elements, 2*node_dof)
    back_dof_lin = collect(node_dof-1:-1:0)
    for i=1:n_elements
        uid = T[i,1]
        vid = T[i,2]
        id_map[i, 1:node_dof] = uid * node_dof * ones(node_dof) - back_dof_lin
        id_map[i, node_dof+1:2*node_dof] = vid * node_dof * ones(node_dof) - back_dof_lin
    end

    # truss
    ex_id = [1, 4]
    xy_id = [1, 2, 4, 5]

    I = Int[]
    J = Int[]
    V = Float64[]
    # K ∈ 2 n_node x 2 n_node
    for e=1:n_elements
        u = X[T[e,1],:]
        v = X[T[e,2],:]
        L = norm(u - v)

        # element local stiffness matrix
        K_le = local_stiffness_matrix(L, E, A[e])

        # local -> global rotation matrix
        R_b = local_coord_rotation(u, v)
        R = zeros(Float64, 2*full_node_dof, 2*full_node_dof)
        for k=1:(Int).(full_node_dof/3) * 2
            R[k*3-3+1:k*3, k*3-3+1:k*3] = R_b
        end

        K_Ge = R[ex_id, xy_id]' * K_le * R[ex_id, xy_id]

        for i=1:2*node_dof
            for j=1:2*node_dof
                push!(I, id_map[e,i])
                push!(J, id_map[e,j])
                push!(V, K_Ge[i,j])
            end
        end
    end
    return sparse(I, J, V, sys_dof, sys_dof)
end

function assemble_load_vector(Load::Matrix{Float64}, n_nodes::Int, node_dof::Int)
    @assert(size(Load,2) == node_dof + 1)
    sys_dof::Int = node_dof * n_nodes
    n_load_nodes = size(Load, 1)
    I = Int[]
    V = Float64[]
    for i=1:n_load_nodes
        v_id = Int.(Load[i,1])
        append!(I, node_dof*(v_id-1)+1 : node_dof*v_id)
        append!(V, Load[i, 2:1:node_dof+1])
    end
    return sparsevec(I, V, sys_dof)
end

function dof_permutation(S::Matrix{Int}, n_nodes::Int, node_dof::Int)
    sys_dof::Int = node_dof * n_nodes
    n_fix = size(S,1)
    n_dof_fix = sum(sum(S[:,2:1+node_dof] .== 1))

    full_dof_fix = zeros(Bool, sys_dof)
    for i=1:n_fix
        v_id = Int.(S[i,1])
        full_dof_fix[node_dof*(v_id-1)+1 : node_dof*v_id] = S[i, 2:1+node_dof]
        #TODO: node_dof here might be incompatible
    end

    perm = zeros(Int, sys_dof)
    n_dof_free = sys_dof - n_dof_fix
    free_tail = 1
    fix_tail = n_dof_free + 1
    for i=1:sys_dof
        if !full_dof_fix[i]
            perm[free_tail] = i
            free_tail += 1
        else
            perm[fix_tail] = i
            fix_tail += 1
        end
    end
    I = collect(1:sys_dof)
    perm_spm = sparse(I, perm, ones(sys_dof), sys_dof, sys_dof)
    return perm, perm_spm, n_dof_free
end

function weight_calculation()
    # assemble stiffness matrix

    # solve linear system

    # TODO: buckling sizing is ignored for now
    # calculate ∑ F_e * l_e

    # return W
end
