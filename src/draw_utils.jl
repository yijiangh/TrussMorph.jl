using Makie

function draw_truss!(scene, X::Matrix{Float64}, T::Matrix{Int}, S::Matrix{Int}, area::Vector{Float64}; line_width::Float64=1.0, draw_supp::Bool=true, supp_scale::Float64=0.1, color=:black, xaxis_label::String="x")

    @assert(size(T, 1) == length(area))
    a = reshape([area area]', 2*length(area))
    # a ./= maximum(a)
    a .*= line_width

    # seg_ids = reshape([T[:,1] T[:,2]], length(T[:,1])+length(T[:,2]))
    seg_ids = reshape(T', prod(size(T)))
    max_xlim = maximum(X[:,1]) - minimum(X[:,1])
    max_ylim = maximum(X[:,2]) - minimum(X[:,2])
    max_lim = max(max_xlim, max_ylim)
    limits = FRect(minimum(X[:,1]), minimum(X[:,2]),
                   max_lim, max_lim)

    linesegments!(scene, X[seg_ids, 1], X[seg_ids, 2], linewidth = a,
                  color = color,
                  limits = limits,
                  axis = (names = (axisnames = (xaxis_label, "y"),),
                          grid = (linewidth = (1, 1),),
                          )
                  )
    if draw_supp
        fix_ids = (Int).(S[:,1])
        # arrows!(scene,
        #       vcat(truss.X[fix_ids, 1],truss.X[fix_ids, 1]),
        #       vcat(truss.X[fix_ids, 2],truss.X[fix_ids, 2]),
        #       vcat(truss.S[:,2], zeros(size(truss.S,1))).*supp_scale,
        #       vcat(zeros(size(truss.S,1)), truss.S[:,3]).*supp_scale,
        #       linecolor = :green, arrowcolor = :green)
        scatter!(scene, X[fix_ids, 1], X[fix_ids, 2], color = :green, limits = limits)
    end
    return scene
end

function draw_load!(scene, X::Matrix{Float64}, load::Matrix{Float64}; load_scale::Float64=0.1, xaxis_label::String="x")
    fix_ids = (Int).(load[:,1])
    max_xlim = maximum(X[:,1]) - minimum(X[:,1])
    max_ylim = maximum(X[:,2]) - minimum(X[:,2])
    max_lim = max(max_xlim, max_ylim)
    limits = FRect(minimum(X[:,1]), minimum(X[:,2]),
                   max_lim, max_lim)
    arrows!(scene, X[fix_ids, 1], X[fix_ids, 2],
            load[:,2].*load_scale, load[:,3].*load_scale,
            linecolor = :pink, arrowcolor = :pink, limits = limits,
            axis = (names = (axisnames = (xaxis_label, "y"),),)
            )
end

function draw_deformed!(scene, truss::Truss, U::Vector{Float64}, node_dof::Int;
    line_width::Float64=10, deform_scale::Float64=5.0)

    n_v = size(truss.X, 1)
    n_e = size(truss.T, 1)
    a = ones(Float64, 2*n_e)
    a .*= line_width

    seg_ids = reshape([truss.T[:,1] truss.T[:,2]], length(truss.T[:,1])+length(truss.T[:,2]))
    max_lim = max(maximum(truss.X[:,1]) - minimum(truss.X[:,1]),
                  maximum(truss.X[:,2]) - minimum(truss.X[:,2]))
    limits = FRect(minimum(truss.X[:,1]), minimum(truss.X[:,2]),
                   max_lim, max_lim)

    node_U = reshape(U, (node_dof, Int.(length(U)/node_dof)))'
    DX = truss.X[seg_ids, 1] + node_U[seg_ids, 1]
    DY = truss.X[seg_ids, 2] + node_U[seg_ids, 2]
    linesegments!(scene, DX .* deform_scale, DY .* deform_scale,
                  linewidth = a,
                  color = :blue,
                  limits = limits,
                  # axis = (names = (axisnames = ("x", "y"),),
                  #         grid = (linewidth = (1, 1),),
                  #         )
                          )
end
