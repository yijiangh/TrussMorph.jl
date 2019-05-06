using Makie

function draw_truss!(scene, truss::Truss, area::Vector{Float64}; line_width=10, draw_supp=true, supp_scale=0.1)
    @assert(size(truss.T, 1) == length(area))
    a = reshape([area area]', 2*length(area))
    a ./= maximum(a)
    a .*= line_width

    seg_ids = reshape([truss.T[:,1] truss.T[:,2]], length(truss.T[:,1])+length(truss.T[:,2]))
    max_lim = max(maximum(truss.X[:,1]) - minimum(truss.X[:,1]),
                  maximum(truss.X[:,2]) - minimum(truss.X[:,2]))
    limits = FRect(minimum(truss.X[:,1]), minimum(truss.X[:,2]),
                   max_lim, max_lim)

    linesegments!(scene, truss.X[seg_ids, 1], truss.X[seg_ids, 2], linewidth = a,
                  color = :black,
                  limits = limits,
                  axis = (names = (axisnames = ("x", "y"),),
                          grid = (linewidth = (1, 1),),
                          ))
    if draw_supp
        fix_ids = (Int).(truss.S[:,1])
        # arrows!(scene,
        #       vcat(truss.X[fix_ids, 1],truss.X[fix_ids, 1]),
        #       vcat(truss.X[fix_ids, 2],truss.X[fix_ids, 2]),
        #       vcat(truss.S[:,2], zeros(size(truss.S,1))).*supp_scale,
        #       vcat(zeros(size(truss.S,1)), truss.S[:,3]).*supp_scale,
        #       linecolor = :green, arrowcolor = :green)
        scatter!(scene, truss.X[fix_ids, 1], truss.X[fix_ids, 2], color = :green, limits = limits)
    end
end

function draw_load!(scene, truss::Truss, load::Matrix{Float64}; load_scale::Float64=0.1)
    fix_ids = (Int).(load[:,1])
    max_lim = max(maximum(truss.X[:,1]) - minimum(truss.X[:,1]),
                  maximum(truss.X[:,2]) - minimum(truss.X[:,2]))
    limits = FRect(minimum(truss.X[:,1]), minimum(truss.X[:,2]),
                   max_lim, max_lim)
    arrows!(scene, truss.X[fix_ids, 1], truss.X[fix_ids, 2],
            load[:,2].*load_scale, load[:,3].*load_scale,
            linecolor = :pink, arrowcolor = :pink, limits = limits)
end
