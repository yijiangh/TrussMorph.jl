using TrussMorph
tm = TrussMorph
using Test
using Makie
using Colors

dir = "/Users/yijiangh/Dropbox (MIT)/Course_Work/6.838_2019_spring/final_project/code/gh_validation"
# fp0 = joinpath(dir, "wiggle_truss_min_w.json")
# fp1 = joinpath(dir, "wiggle_truss_min_d.json")
# load_fp = joinpath(dir, "wiggle_truss_min_d_load_case.json")

fp0 = joinpath(dir, "2D_truss.json")
fp1 = joinpath(dir, "2D_truss_1.json")
load_fp = joinpath(dir, "2D_truss_load_case.json")

node_dof = 3 # 2 for truss model, 3 for frames
full_node_dof = 3 # fixed for 2D cases

t0 = parse_truss_json(fp0)
t1 = parse_truss_json(fp1)
load = parse_load_json(load_fp, node_dof)

design_var_ids = [4, 6]
# design_var_ids = [3, 4, 5, 6]

# full_design_var_ids = vcat(collect(3:12), collect(15:24))
# full_design_var_ids = reshape(full_design_var_ids, (2, Int.(length(full_design_var_ids) / 2)))'
# design_var_ids = full_design_var_ids[:,2]

path_disc = 10
parm_weight = 1.0
parm_smooth = 100.0

morph_path, ptwise_smoothness_energy, ptwise_weight_energy, ptwise_total_energy, parm =
    tm.compute_morph_path(t0, t1, load, node_dof, full_node_dof, design_var_ids, path_disc=path_disc,parm_smooth=parm_smooth, parm_weight=parm_weight)

scene = Scene()
tm.draw_load!(scene, t0, load, load_scale=0.002)

plen = length(morph_path)
r = ones(Float64, size(t0.T,1))
color_array = Array{RGBAf0}(undef, plen)

for i=1:plen
    mcolor = Float32.((plen - i) / plen) * RGBAf0(1.0,0.0,0.0,0.1) + Float32.(i / plen) * RGBAf0(0.0,0.0,1.0,1)
    color_array[i] = mcolor

    tm.draw_truss!(scene, morph_path[i], t0.T, t0.S, r, supp_scale=0.2, color=mcolor)
    # tm.draw_truss!(scene, t0.X, t0.T, t0.S, r, supp_scale=0.2, color=mcolor)

    display(scene)
    sleep(0.5)
end

# TODO: plot weight evolution
scene_sE = barplot(ptwise_smoothness_energy, color=color_array,
                   axis = (names = (axisnames = ("iter", "smoothness (neighbor) energy"),),
                   grid = (linewidth = (1, 1),),
                   ))
scene_wE = barplot(ptwise_weight_energy, color=color_array,
                  axis = (names = (axisnames = ("iter", "weight energy"),),
                  grid = (linewidth = (1, 1),),
                  ))
scene_totalE = barplot(ptwise_total_energy, color=color_array,
                 axis = (
                    names = (axisnames = ("iter", string("total energy- pS:", parm[1], ", pW:", parm[2])),),
                    grid = (linewidth = (1, 1),),
                 )
                 )

scene_final = hbox(
    vbox(scene_sE, scene_wE, scene_totalE),
    scene,
    )

display(scene_final)

# tm.draw_deformed!(scene, t, U, node_dof)
