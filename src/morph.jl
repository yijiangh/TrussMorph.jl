using TrussMorph
tm = TrussMorph
using Optim
using Printf

function compute_morph_path(t0::tm.Truss, t1::tm.Truss, load::Matrix{Float64},
    node_dof::Int, full_node_dof::Int, design_var_ids::Vector{Int}; path_disc::Int=5, parm_smooth::Float64=100.0, parm_weight::Float64=50.0, do_opt::Bool=true, kwargs...)

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

    var_chuck = length(design_var_ids)::Int
    X0_var = reshape(t0.X', prod(size(t0.X)))[design_var_ids]
    X1_var = reshape(t1.X', prod(size(t1.X)))[design_var_ids]

    # TODO: symmetric variable
    weight_fn = tm.get_weight_calculation_fn(t0.X, r, t0.T, F_m, perm, n_dof_free,
    node_dof, full_node_dof, t0.mp, design_var_ids)

    # set up Xpath0
    Xpath0 = Vector{Float64}(undef, path_disc * var_chuck)
    for i=1:path_disc
        Xpath0[(i-1)*var_chuck + 1 : i*var_chuck] =
            Float64.((path_disc - i + 1) / (path_disc + 1)) * X0_var +
            Float64.(i / (path_disc + 1)) * X1_var
    end

    L_I = collect(1:path_disc+1)
    L_J = collect(2:path_disc+2)
    L_V = ones(Float64, path_disc+1)
    L = sparse(L_I, L_J, L_V, path_disc+2, path_disc+2)
    L = L + L'
    rowsums = vec(sum(L,dims=2))
    L = L - spdiagm(0 => rowsums)

    function path_energy(Xpath::Vector{Float64})
        Xmat = vcat(X0_var', reshape(Xpath, (var_chuck, path_disc))', X1_var')
        dXpath_dt = Xmat[2:end, :] - Xmat[1:end-1, :]
        path_weight = 0
        for i=1:path_disc
            path_weight += weight_fn(Xmat[1+i,:])
        end
        path_weight_vec = zeros(Float64, size(Xmat, 1))
        for i=1:size(Xmat, 1)
            path_weight_vec[i] = weight_fn(Xmat[i,:])
        end
        # ddt = (L * path_weight_vec).^2

        # ∑ |dϕ/dt|^2 + |dWF/dt|^2
        return parm_smooth * sum(dXpath_dt.^2) + parm_weight * sum((path_weight_vec[2:end] - path_weight_vec[1:end-1]).^2)

        # ∑ |dWF/dt|^2
        # return sum((path_weight_vec[2:end] - path_weight_vec[1:end-1]).^2)

        # ∑ |d^2WF/dt^2|^2
        # return sum(ddt[2:end-1])
        # return sum(ddt)

        # ∑ |dϕ/dt|^2 + |WF(X(t))|
        # return parm_smooth * sum(dXpath_dt.^2) + parm_weight * path_weight
    end

    @show init_sumE = path_energy(Xpath0)
    opt_sumE = init_sumE

    # run opt
    if do_opt
        if :box_constraint in keys(kwargs) && kwargs[:box_constraint] != undef
            bc_var = kwargs[:box_constraint]
            @assert(size(bc_var, 1) == length(design_var_ids))
            @assert(size(bc_var, 2) == 2)
            println("using box constraint...")

            Xpath_bc_lower = reshape(bc_var[:,1] * ones(path_disc)', path_disc * var_chuck)
            Xpath_bc_upper = reshape(bc_var[:,2] * ones(path_disc)', path_disc * var_chuck)
            # @assert(length(Xpath_bc_lower) == length(Xpath0))
            # @assert(length(Xpath_bc_upper) == length(Xpath0))

            # inner_optimizer = GradientDescent()
            # res = Optim.optimize(path_energy, box_constraint[:,1], box_constraint[:,2], Xpath0, Fminbox(inner_optimizer))
            res = Optim.optimize(path_energy, Xpath_bc_lower, Xpath_bc_upper, Xpath0)
        else
            println("unconstrained...")
            # Optim constructs an approximate gradient for us using central finite differencing
            res = Optim.optimize(path_energy, Xpath0, LBFGS())
        end

        @show Optim.summary(res)
        @show Optim.iterations(res)
        @show Optim.iteration_limit_reached(res)
        @show Optim.f_calls(res)
        @show Optim.converged(res)

        X_var_star_vec = Optim.minimizer(res)
        opt_sumE = Optim.minimum(res)
        @show opt_sumE
    else
        X_var_star_vec = Xpath0
    end

    X_var_star = reshape(X_var_star_vec, (var_chuck, path_disc))'
    X_var_init = reshape(Xpath0, (var_chuck, path_disc))'

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

    for i=2:path_disc+1
        morph_path[i] = copy(X_template)
        for j=1:var_chuck
            morph_path[i][design_var_id_map[j,1], design_var_id_map[j,2]] = X_var_star[i-1, j]
        end
    end

    init_morph_path = Array{Matrix{Float64}}(undef, path_disc + 2)
    init_morph_path[1] = t0.X
    init_morph_path[end] = t1.X

    for i=2:path_disc+1
        init_morph_path[i] = copy(X_template)
        for j=1:var_chuck
            init_morph_path[i][design_var_id_map[j,1], design_var_id_map[j,2]] = X_var_init[i-1, j]
        end
    end

    # function map_full_X_to_design_vars(origX::Matrix{Float64})::Vector{Float64}
    #     return reshape(origX', prod(size(origX)))[design_var_ids]
    # end
    function cal_ptwise_energies(X_var_vec::Vector{Float64})
        X_var = reshape(X_var_vec, (var_chuck, path_disc))'
        ptwise_smoothness_energy = Float64[]
        ptwise_weight_energy = Float64[]
        ptwise_total_energy = Float64[]
        push!(ptwise_smoothness_energy, 0.5 * sum((X0_var - X_var[1,:]).^2))
        push!(ptwise_weight_energy, weight_fn(X0_var))
        push!(ptwise_total_energy, parm_smooth * ptwise_smoothness_energy[end] +
                                   parm_weight * ptwise_weight_energy[end])

        for i=1:path_disc
            push!(ptwise_weight_energy, weight_fn(X_var[i, :]))
            if 1 == i
                prevX_var = X0_var
            else
                prevX_var = X_var[i-1,:]
            end
            if path_disc == i
                nextX_var = X1_var
            else
                nextX_var = X_var[i+1,:]
            end
            push!(ptwise_smoothness_energy, 0.5 * (sum((X_var[i, :] - prevX_var).^2) + sum((X_var[i, :] - nextX_var).^2)))
            push!(ptwise_total_energy, parm_smooth * ptwise_smoothness_energy[end] + parm_weight * ptwise_weight_energy[end])
        end

        push!(ptwise_smoothness_energy, 0.5 * sum((X1_var - X_var[end,:]).^2))
        push!(ptwise_weight_energy, weight_fn(X1_var))
        push!(ptwise_total_energy, parm_smooth * ptwise_smoothness_energy[end] +
                                   parm_weight * ptwise_weight_energy[end])
        return ptwise_smoothness_energy, ptwise_weight_energy, ptwise_total_energy
    end

    opt_sE, opt_wE, opt_totE = cal_ptwise_energies(X_var_star_vec)
    init_sE, init_wE, init_totE = cal_ptwise_energies(Xpath0)

    # return intial morph path and energy as a comparision
    return morph_path, init_morph_path, opt_sE, opt_wE, opt_totE, init_sE, init_wE, init_totE, opt_sumE, init_sumE, [parm_smooth, parm_weight]
end
