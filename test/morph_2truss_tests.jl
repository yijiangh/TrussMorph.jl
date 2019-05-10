using TrussMorph
tm = TrussMorph
using Test
using Makie
using Colors

node_dof = 3 # 2 for truss model, 3 for frames
full_node_dof = 3 # fixed for 2D cases

dir = "/Users/yijiangh/Dropbox (MIT)/Course_Work/6.838_2019_spring/final_project/code/gh_validation"

fp0 = joinpath(dir, "wiggle_truss_min_w.json")
fp1 = joinpath(dir, "wiggle_truss_min_d.json")
load_fp = joinpath(dir, "wiggle_truss_min_d_load_case.json")
full_design_var_ids = vcat(collect(3:12), collect(15:24))
full_design_var_ids = reshape(full_design_var_ids, (2, Int.(length(full_design_var_ids) / 2)))'
design_var_ids = full_design_var_ids[:,2]

# fp0 = joinpath(dir, "2D_truss.json")
# fp1 = joinpath(dir, "2D_truss_1.json")
# load_fp = joinpath(dir, "2D_truss_load_case.json")
# design_var_ids = [4, 6]
# design_var_ids = [3, 4, 5, 6]

t0 = parse_truss_json(fp0)
t1 = parse_truss_json(fp1)
load = parse_load_json(load_fp, node_dof)

path_disc = 15
parm_weight = 10.0
parm_smooth = 100.0

morph_path, init_morph_path, opt_sE, opt_wE, opt_totE, init_sE, init_wE, init_totE, opt_sumE, init_sumE, parm =
    tm.compute_morph_path(t0, t1, load, node_dof, full_node_dof, design_var_ids, path_disc=path_disc,parm_smooth=parm_smooth, parm_weight=parm_weight)

scene = Scene()
init_scene= Scene()
tm.draw_load!(scene, t0, load, load_scale=0.002, xaxis_label="x - OPT")
tm.draw_load!(init_scene, t0, load, load_scale=0.002, xaxis_label="x - INIT(linear)")

plen = length(morph_path)
r = ones(Float64, size(t0.T,1)) * 0.5
color_array = Array{RGBAf0}(undef, plen)

for i=1:plen
    mcolor = Float32.((plen - i) / plen) * RGBAf0(1.0,0.0,0.0,0.1) + Float32.(i / plen) * RGBAf0(0.0,0.0,1.0,1)
    color_array[i] = mcolor

    tm.draw_truss!(scene, morph_path[i], t0.T, t0.S, r, supp_scale=0.2, color=mcolor, xaxis_label="x - OPT")
    tm.draw_truss!(init_scene, init_morph_path[i], t0.T, t0.S, r, supp_scale=0.2, color=mcolor, xaxis_label="x - INIT(linear)")

    display(scene)
    sleep(0.5)
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
                   axis = (names = (axisnames = ("iter INIT", "smoothness (neighbor) energy"),),
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
                    names = (axisnames = (string("iter INIT - total E: ", init_sumE), string("total energy- pS:", parm[1], ", pW:", parm[2])),),
                        ),
                  limits=totE_limits,
                 )

scene_final = hbox(
    vbox(scene_sE, scene_wE, scene_totalE),
    vbox(init_scene_sE, init_scene_wE, init_scene_totalE),
    vbox(scene, init_scene)
    )

display(scene_final)

# tm.draw_deformed!(scene, t, U, node_dof)
