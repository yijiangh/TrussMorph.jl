using Optim

function compute_morph_path(t0::TrussMorph.Truss, t1::TrussMorph.Truss, load::Matrix{Float64},
    node_dof::Int, full_node_dof::Int; path_disc=5)
    @assert(size(t0.X, 1) == size(t1.X, 1))
    @assert(t0.T == t1.T)

    n_v = size(t0.X,1)
    n_e = size(t0.T,1)
    # initial cross secs
    r = ones(size(t0.T,1)) * 0.01
    morph_path = Array{Matrix{Float64}}(undef, path_disc + 2)

    F = tm.assemble_load_vector(load, n_v, node_dof)
    perm_vec, perm, n_dof_free = tm.dof_permutation(t.S, n_v, node_dof)

    F_perm = perm * F
    F_m = Array(F_perm[1:n_dof_free])

    weight = tm.weight_calc(X, r, t.T, F_m, perm, n_dof_free, node_dof, full_node_dof, t.mp)

    return morph_path
end
