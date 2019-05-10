using TrussMorph: SectionProperties, MaterialProperties
using LinearAlgebra: norm, cholesky, I
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
function local_stiffness_matrix(L::Float64, mp::MaterialProperties, sp::SectionProperties)
    E = mp.E
    G = mp.G
    A = sp.A
    I = sp.Iz
    fs = 1.0 # shear deformation const
    beta = (12*E*I*fs)/(G*A*L^2)

    K_beam = zeros(Float64, 4, 4)
    K_beam[1,2] = 6*L;
    K_beam[1,3] = -12;
    K_beam[1,4] = 6*L;
    K_beam[2,3] = -6*L;
    K_beam[2,4] = L^2*(2-beta);
    K_beam[3,4] = -6*L;
    K_beam = K_beam' + K_beam;

    K_beam += Diagonal([12, L^2*(4+beta), 12, L^2*(4+beta)])
    K_beam *= (E*I)/(L^3*(1+beta))

    K_le = zeros(Float64, 6, 6)
    K_le[1,1] = E * A / L
    K_le[1,4] = -E * A / L
    K_le[4,1] = -E * A / L
    K_le[4,4] = E * A / L
    K_le[2:3,2:3] = K_beam[1:2,1:2];
    K_le[2:3,5:6] = K_beam[1:2,3:4];
    K_le[5:6,2:3] = K_beam[3:4,1:2];
    K_le[5:6,5:6] = K_beam[3:4,3:4];
    return K_le
end

# TODO: can be optimized using symbolic substitution, only X is changing here
# NOTE: r might be set to be uniform, if we are using ∑ f_e * l_e formula?
function assemble_global_stiffness_matrix(X::Matrix{Float64}, T::Matrix{Int64},
                                          r::Vector{Float64}, mp::MaterialProperties,
                                          node_dof::Int, full_node_dof::Int)
    n_elements = size(T,1)
    n_nodes = size(X,1)

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

    # cross sec properties
    sp_vec = Array{SectionProperties}(undef, n_elements)
    for e=1:n_elements
        sp_vec[e] = compute_round_section_properties(r[e])
    end

    # for elemental force calc
    KR_es = Array{Matrix{Float64},1}(undef, n_elements)

    if 2 == node_dof
        ex_id = [1, 4]
        xy_id = [1, 2, 4, 5]
    else
        ex_id = collect(1:2*node_dof)
        xy_id = collect(1:2*node_dof)
    end

    # TODO: array size node_dof^2 * n_elements
    I = Int[]
    J = Int[]
    V = Float64[]
    # K ∈ (node_dof*n_node)^2
    for e=1:n_elements
        u = X[T[e,1],:]
        v = X[T[e,2],:]
        L = norm(u - v)

        # element local stiffness matrix
        K_le = local_stiffness_matrix(L, mp, sp_vec[e])

        # local -> global rotation matrix
        R_b = local_coord_rotation(u, v)
        R = zeros(Float64, 2*full_node_dof, 2*full_node_dof)
        for k=1:(Int).(full_node_dof/3) * 2
            R[k*3-3+1:k*3, k*3-3+1:k*3] = R_b
        end

        K_Ge = R[ex_id, xy_id]' * K_le[ex_id, ex_id] * R[ex_id, xy_id]
        KR_es[e] = K_le * R[ex_id, xy_id]

        for i=1:2*node_dof
            for j=1:2*node_dof
                push!(I, id_map[e,i])
                push!(J, id_map[e,j])
                push!(V, K_Ge[i,j])
            end
        end
    end
    return sparse(I, J, V, sys_dof, sys_dof), KR_es, id_map
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
        #TODO: 3D node_dof here might be incompatible, 2D is okay
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

function get_weight_calculation_fn(X_full::Matrix{Float64}, r::Vector{Float64}, T::Matrix{Int64},
    F_perm_m::Vector{Float64}, perm::SparseMatrixCSC{Float64, Int}, n_dof_free::Int, node_dof::Int, full_node_dof::Int, mp::MaterialProperties, design_var_id::Vector{Int}; opt_compliance::Bool=false)::Function

    # assemble stiffness matrix
    n_v = size(X_full, 1)
    n_e = size(T, 1)
    sys_dof = n_v * node_dof

    # split design variable + the rest
    X_template = copy(X_full)
    design_var_id_map = zeros(Int, length(design_var_id), 2)
    for i=1:length(design_var_id)
        design_var_id_map[i,1] = Int.((design_var_id[i] + design_var_id[i]%2) / 2)
        design_var_id_map[i,2] = 2 - design_var_id[i]%2
        X_template[design_var_id_map[i,1], design_var_id_map[i,2]] = 0
    end
    # @show design_var_id
    # @show X_template
    # @show design_var_id_map

    function calc_weight(X_var::Vector{Float64})
        X_new = copy(X_template) #this is bad
        for i=1:length(design_var_id)
            X_new[design_var_id_map[i,1], design_var_id_map[i,2]] = X_var[i]
        end

        K, KR_es, id_map = assemble_global_stiffness_matrix(X_new, T, r, mp, node_dof, full_node_dof)
        K_perm = perm * K * perm'
        K_mm = K_perm[1:n_dof_free, 1:n_dof_free]

        # solve linear system
        # U_m = K_mm\F_perm_m
        K_mm = (K_mm + K_mm')/2
        Kpf = cholesky(K_mm, check=false)
        if !LinearAlgebra.issuccess(Kpf)
            return 1e10
        end
        U_m = Kpf\F_perm_m

        # U_I, U_V = findnz(U_m)
        U_perm = vcat(U_m, zeros(sys_dof - n_dof_free))

        # compliance
        if opt_compliance
            return F_perm_m' * U_m
        end

        U = perm' * U_perm

        # TODO: buckling sizing is ignored for now
        # calculate ∑ F_e * l_e
        weight = 0.0
        if node_dof == 2
            e_react_dof = 1 # truss
        else
            e_react_dof = 3 # frame
        end

        eF = zeros(n_e, e_react_dof*2)
        for e=1:n_e
            eF[e,:] = KR_es[e] * U[id_map[e,:]]
            eL = norm(X_new[T[e,1], :] - X_new[T[e,2], :])
            weight += abs(eF[e,:][1]) * eL
        end
        # @show eF
        # @show U
        return weight
    end

    return calc_weight
end
