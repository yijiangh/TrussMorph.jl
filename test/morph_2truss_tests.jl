using TrussMorph
tm = TrussMorph
using Test
using Makie
using Colors
using Dates
using Printf

recompute = true
plot = true
node_dof = 3 # 2 for truss model, 3 for frames
full_node_dof = 3 # fixed for 2D cases
box_constr = undef

dir = "/Users/yijiangh/Dropbox (MIT)/Course_Work/6.838_2019_spring/final_project/code/gh_validation"
result_file_dir = joinpath(pwd(),"test","results")

# st_file_name  = "wiggle_truss_min_d.json"
# end_file_name  = "wiggle_truss_min_w.json"
# fp0 = joinpath(dir, st_file_name)
# fp1 = joinpath(dir, end_file_name)
# load_fp = joinpath(dir, "wiggle_truss_min_d_load_case.json")
# full_design_var_ids = vcat(collect(3:12), collect(15:24))
# full_design_var_ids = reshape(full_design_var_ids, (2, Int.(length(full_design_var_ids) / 2)))'
# design_var_ids = full_design_var_ids[:,2]

# st_file_name  = "funicular_arch.json"
# end_file_name  = "funicular_cable.json"
# fp0 = joinpath(dir, st_file_name)
# fp1 = joinpath(dir, end_file_name)
# load_fp = joinpath(dir, "funicular_arch_load_case.json")
# design_var_ids = collect(4:2:20)

st_file_name  = "2D_truss_0.json"
end_file_name  = "2D_truss_1.json"
fp0 = joinpath(dir, st_file_name)
fp1 = joinpath(dir, end_file_name)
load_fp = joinpath(dir, "2D_truss_0_load_case.json")
design_var_ids = [4, 6]
# design_var_ids = [3, 4, 5, 6]

# st_file_name  = "cm_truss_0.json"
# end_file_name  = "cm_truss_1.json"
# fp0 = joinpath(dir, st_file_name)
# fp1 = joinpath(dir, end_file_name)
# load_fp = joinpath(dir, "cm_truss_0_load_case.json")
# # design_var_ids = [4, 6, 8]
# design_var_ids = [3, 4, 6, 7, 8]
# X0_var = reshape(t0.X', prod(size(t0.X)))[design_var_ids]
# X1_var = reshape(t1.X', prod(size(t1.X)))[design_var_ids]
# box_constr = [0.2 3; -5 5; -5 5; 3 5.7; -5 5]
# box_constr += hcat(X0_var, X0_var)

path_disc = 20
parm_weight = 1.0
parm_smooth = 1e7

# plot parameters
load_scale = 0.1
line_width = 4.0
supp_scale = 0.3

t0,_ = parse_truss_json(fp0)
t1,_ = parse_truss_json(fp1)
load = parse_load_json(load_fp, node_dof)

morph_path, init_morph_path, opt_sE, opt_wE, opt_totE, init_sE, init_wE, init_totE, opt_sumE, init_sumE, parm =
    tm.compute_morph_path(t0, t1, load, node_dof, full_node_dof, design_var_ids, path_disc=path_disc,parm_smooth=parm_smooth, parm_weight=parm_weight, do_opt=recompute, box_constraint=box_constr)

pure_st_file_name = SubString(st_file_name, 1:length(st_file_name)-length(".json"))
pure_end_file_name = SubString(end_file_name, 1:length(end_file_name)-length(".json"))
result_file_name = pure_st_file_name * "-" *  pure_end_file_name
f_file_dir = joinpath(result_file_dir, result_file_name)

if !recompute
    searchdir(path,key) = filter(x->occursin(key,x), readdir(path))

    opt_sE = Float64[]
    opt_wE = Float64[]
    opt_totE = Float64[]

    @assert(ispath(f_file_dir))
    result_jsons = searchdir(f_file_dir,string("-", path_disc+2, ".json"))
    # @show result_jsons
    morph_path = Matrix{Float64}[]
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
else
    tm.save_morph_path_json(morph_path, result_file_dir, st_file_name, end_file_name, t0, opt_sE, opt_wE, opt_totE, opt_sumE, parm_smooth, parm_weight)
end

if plot
    scene = Scene()
    init_scene= Scene()
    tm.draw_load!(scene, t0.X, load, load_scale=load_scale, xaxis_label="x - OPT")
    tm.draw_load!(init_scene, t0.X, load, load_scale=load_scale, xaxis_label="x - INIT(linear)")

    plen = length(morph_path)
    r = ones(Float64, size(t0.T,1)) * 0.5
    color_array = Array{RGBAf0}(undef, plen)

    for i=1:plen
        mcolor = Float32.((plen - i) / plen) * RGBAf0(1.0,0.0,0.0,0.1) + Float32.(i / plen) * RGBAf0(0.0,0.0,1.0,1)
        color_array[i] = mcolor

        tm.draw_truss!(scene, morph_path[i], t0.T, t0.S, r, supp_scale=supp_scale, color=mcolor, xaxis_label="x - OPT", line_width=line_width)

        tm.draw_truss!(init_scene, init_morph_path[i], t0.T, t0.S, r, supp_scale=supp_scale, color=mcolor, xaxis_label="x - INIT(linear)", line_width=line_width)
        # display(scene)
        # sleep(0.2)
    end

    sE_limits = FRect(1, 0, plen, max(maximum(opt_sE), maximum(init_sE)))
    wE_limits = FRect(1, 0, plen, max(maximum(opt_wE), maximum(init_wE)))
    totE_limits = FRect(1, 0, plen, max(maximum(opt_totE), maximum(init_totE)))

    scene_sE = barplot(opt_sE, color=color_array,
                       axis = (names = (axisnames = ("iter - OPT", "smoothness (neighbor) energy"),),
                              ),
                       limits=sE_limits,
                       )
    scene_wE = barplot(opt_wE, color=color_array,
                      axis = (names = (axisnames = ("iter - OPT", "weight energy"),),
                             ),
                      limits=wE_limits,
                      )
    scene_totalE = barplot(opt_totE, color=color_array,
                     axis = (
                        names = (axisnames = (string("iter - OPT total E: ", opt_sumE), string("total energy- pS:", parm[1], ", pW:", parm[2])),),
                            ),
                     limits=totE_limits,
                     )

    init_scene_sE = barplot(init_sE, color=color_array,
                       axis = (
                            names = (axisnames = ("iter INIT", "smoothness (neighbor) energy"),),
                              ),
                       limits=sE_limits,
                       )
    init_scene_wE = barplot(init_wE, color=color_array,
                      axis = (names = (axisnames = ("iter INIT", "weight energy"),),
                             ),
                      limits=wE_limits,
                      )
    init_scene_totalE = barplot(init_totE, color=color_array,
                     axis = (
                        names = (axisnames = (string("iter INIT -total E: ", init_sumE), string("total energy- pS:", parm[1], ", pW:", parm[2])),),
                            ),
                      limits=totE_limits,
                     )

    scene_final = hbox(
        vbox(scene_sE, scene_wE, scene_totalE),
        vbox(init_scene_sE, init_scene_wE, init_scene_totalE),
        vbox(scene, init_scene)
        )

    Makie.save(joinpath(f_file_dir, result_file_name * "_" * string(Dates.now()) * ".png"), scene_final)

    anim_sc = Scene()
    # save morph gif
    # record(anim_sc, joinpath(f_file_dir, result_file_name * "_" * string(Dates.now()) * ".gif"), 1:plen) do i
    #     tm.draw_truss!(anim_sc, morph_path[i], t0.T, t0.S, r, supp_scale=supp_scale, color=color_array[i], xaxis_label="x - OPT", line_width)
    #     display(anim_sc)
    #     sleep(1)
    # end
    # t = Node(i)
    # color=color_array[lift(i->i,t)]

    record(anim_sc, joinpath(f_file_dir, result_file_name * "_" * string(Dates.now()) * ".gif"), 1:plen) do i
        anim_sc = Scene()
        mcolor = Float32.((plen - i) / plen) * RGBAf0(1.0,0.0,0.0,1) + Float32.(i / plen) * RGBAf0(0.0,0.0,1.0,1)

        tm.draw_load!(anim_sc, morph_path[i], load, load_scale=load_scale, xaxis_label="x - OPT")
        tm.draw_truss!(anim_sc, morph_path[i], t0.T, t0.S, r, supp_scale=supp_scale, xaxis_label="x - OPT", color = mcolor, line_width=line_width)

        result_img_name = result_file_name * string("_", @sprintf("%03d",i),"-",plen) * ".png"
        Makie.save(joinpath(f_file_dir, result_img_name), anim_sc)

        display(anim_sc)
        sleep(0.1)
    end

    display(scene_final)

end
