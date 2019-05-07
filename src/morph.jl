using Optim
# include("common_types.jl")

function compute_morph_path(t0::TrussMorph.Truss, t1::TrussMorph.Truss, load::Matrix{Float64},
    node_dof::Int, full_node_dof::Int; path_disc=5)
    @assert(size(t0.X, 1) == size(t1.X, 1))
    @assert(t0.T == t1.T)

    n_v = size(t0.X,1)
    n_e = size(t0.T,1)
    # initial cross secs
    A = ones(size(t0.T,1))* pi * 0.01^2
    sys_dof = node_dof * n_v

    morph_path = Array{Matrix{Float64}}(undef, path_disc + 2)

    return morph_path
end
