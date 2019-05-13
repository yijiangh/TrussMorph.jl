using TrussMorph
tm = TrussMorph
using Test
using Makie
using Colors
# using ForwardDiff
# using ReverseDiff

parse_result = true

node_dof = 3 # 2 for truss model, 3 for frames
full_node_dof = 3 # fixed for 2D cases

dir = "/Users/yijiangh/Dropbox (MIT)/Course_Work/6.838_2019_spring/final_project/code/gh_validation"
result_file_dir = joinpath(pwd(),"test","results")

st_file_name  = "2D_truss.json"
end_file_name  = "2D_truss_1.json"
fp0 = joinpath(dir, st_file_name)
fp1 = joinpath(dir, end_file_name)
load_fp = joinpath(dir, "2D_truss_load_case.json")
design_var_ids = [4, 6]

path_disc = 20

t0,_ = parse_truss_json(fp0)
t1,_ = parse_truss_json(fp1)
load = parse_load_json(load_fp, node_dof)

n_v = size(t0.X,1)
n_e = size(t0.T,1)

# initial cross secs
r = ones(size(t0.T,1)) * 0.01
F = tm.assemble_load_vector(load, n_v, node_dof)
perm_vec, perm, n_dof_free = tm.dof_permutation(t0.S, n_v, node_dof)

F_perm = perm * F
F_m = Array(F_perm[1:n_dof_free])

weight_fn = tm.get_weight_calculation_fn(t0.X, r, t0.T, F_m, perm, n_dof_free,
node_dof, full_node_dof, t0.mp, design_var_ids)

if parse_result
    searchdir(path,key) = filter(x->occursin(key,x), readdir(path))

    opt_sE = Float64[]
    opt_wE = Float64[]
    opt_totE = Float64[]
    opt_sumE = 0.0
    morph_path = Matrix{Float64}[]
    parm= zeros(2)

    pure_st_file_name = SubString(st_file_name, 1:length(st_file_name)-length(".json"))
    pure_end_file_name = SubString(end_file_name, 1:length(end_file_name)-length(".json"))
    result_file_name = pure_st_file_name * "-" *  pure_end_file_name
    f_file_dir = joinpath(result_file_dir, result_file_name)
    if ispath(f_file_dir)
        result_jsons = searchdir(f_file_dir,string("-", path_disc+2, ".json"))
        # @show result_jsons
        morph_data = Dict()
        for json_fp in result_jsons
            full_json_fp = joinpath(f_file_dir, json_fp)

            parsed_truss, morph_data = parse_truss_json(full_json_fp, parse_morph=true)
            push!(morph_path, parsed_truss.X)
            push!(opt_wE, morph_data["weight_energy"])
            push!(opt_sE, morph_data["smoothness_energy"])
            push!(opt_totE, morph_data["tot_energy"])
            parm = [morph_data["smooth_parameter"], morph_data["weight_parameter"]]
            opt_sumE = morph_data["sum_energy"]
        end
        # @show opt_sumE
        # @show parsed_morph_path
    end

    # X_mat -> x_var
    design_var_id_map = zeros(Int, length(design_var_ids), 2)
    for i=1:length(design_var_ids)
        design_var_id_map[i,1] = Int.((design_var_ids[i] + design_var_ids[i]%2) / 2)
        design_var_id_map[i,2] = 2 - design_var_ids[i]%2
    end

    morph_var_path = zeros(length(morph_path) , length(design_var_ids))
    for t = 1:length(morph_path)
        for i=1:length(design_var_ids)
            morph_var_path[t,i] = morph_path[t][design_var_id_map[i,1], design_var_id_map[i,2]]
        end
    end
end

sc = Scene()
y1 = -1.2:0.01:1.2
y2 = -1.2:0.01:1.2
wF(x, y) = log(1 + weight_fn([x,y]))
# wF(x, y) = weight_fn([x,y])

p1 = surface!(sc, y1, y2, wF, transparency = parse_result)

# z = Float64[wF(x, y) for x in y1, y in y2]
# wireframe!(sc, y1, y2, z)

if parse_result
    morph_weight = zeros(length(morph_path))
    for i = 1:length(morph_path)
        morph_weight[i] = wF(morph_var_path[i,1], morph_var_path[i,2])
    end

    pts = vec(Point3f0.(morph_var_path[:,1], morph_var_path[:,2], morph_weight))
    scatter!(sc, pts, color=:red, markersize=0.08)
    poly_vec = pts[2:end] - pts[1:end-1]
    push!(poly_vec, [0,0,0])
    arrows!(sc, pts, poly_vec, arrowsize = 0.04, linewidth = 1)
end

# contour3d!(sc, y1, y2, (x, y) -> f(x,y), levels = 15, linewidth = 3)

display(sc)
