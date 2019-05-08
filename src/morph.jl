using TrussMorph
tm = TrussMorph
using Optim

function compute_morph_path(t0::tm.Truss, t1::tm.Truss, load::Matrix{Float64},
    node_dof::Int, full_node_dof::Int; path_disc::Int=5)

    @assert(size(t0.X, 1) == size(t1.X, 1))
    @assert(t0.T == t1.T)

    n_v = size(t0.X,1)
    n_e = size(t0.T,1)

    # initial cross secs
    r = ones(size(t0.T,1)) * 0.01
    F = tm.assemble_load_vector(load, n_v, node_dof)
    perm_vec, perm, n_dof_free = tm.dof_permutation(t0.S, n_v, node_dof)

    F_perm = perm * F
    F_m = Array(F_perm[1:n_dof_free])

    design_var_ids = [4, 6]
    var_chuck = length(design_var_ids)::Int
    X0_var = reshape(t0.X', prod(size(t0.X)))[design_var_ids]
    X1_var = reshape(t1.X', prod(size(t1.X)))[design_var_ids]

    weight_fn = tm.get_weight_calculation_fn(t0.X, r, t0.T, F_m, perm, n_dof_free,
    node_dof, full_node_dof, t0.mp, design_var_ids)

    # set up Xpath0
    Xpath0 = Vector{Float64}(undef, path_disc * var_chuck)
    for i=1:path_disc
        Xpath0[(i-1)*var_chuck + 1 : i*var_chuck] =
            Float64.((path_disc - i + 1) / (path_disc + 1)) * X0_var +
            Float64.(i / (path_disc + 1)) * X1_var
    end

    parm_smooth = 100.0
    parm_weight = 1.0

    function path_energy(Xpath::Vector{Float64})
        Xmat = vcat(X0_var', reshape(Xpath, (var_chuck, path_disc))', X1_var')
        dXpath_dt = Xmat[2:end, :] - Xmat[1:end-1, :]
        path_weight = 0
        for i=1:path_disc
            path_weight += weight_fn(Xmat[1+i,:])
        end
        return parm_smooth * sum(dXpath_dt.^2) + parm_weight * path_weight
    end

    # run opt
    res = Optim.optimize(path_energy, Xpath0, LBFGS())
    @show summary(res)
    X_var_star = Optim.minimizer(res)
    @show Optim.minimum(res)

    # X_var_star = Xpath0
    @show path_energy(Xpath0)

    X_var_star = reshape(X_var_star, (var_chuck, path_disc))'

    X_template = copy(t0.X)
    dim = size(t0.X,2)
    design_var_id_map = zeros(Int, var_chuck, 2)
    for i=1:length(design_var_ids)
        design_var_id_map[i,1] = Int.((design_var_ids[i] + design_var_ids[i]%dim) / dim)
        design_var_id_map[i,2] = dim - design_var_ids[i]%dim
        X_template[design_var_id_map[i,1], design_var_id_map[i,2]] = 0
    end

    morph_path = Array{Matrix{Float64}}(undef, path_disc + 2)
    morph_path[1] = t0.X
    morph_path[end] = t1.X
    # @show morph_path
    # @show design_var_id_map
    # @show X_var_star

    for i=2:path_disc+1
        morph_path[i] = copy(X_template)
        for j=1:var_chuck
            # @show design_var_id_map[j,1]
            # @show design_var_id_map[j,2]
            # @show X_var_star[i-1, j]

            morph_path[i][design_var_id_map[j,1], design_var_id_map[j,2]] = X_var_star[i-1, j]
        end
    end

    return morph_path
end
