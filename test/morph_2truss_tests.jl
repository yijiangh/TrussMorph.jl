using TrussMorph
tm = TrussMorph
using Test
using Makie

fp0 = "/Users/yijiangh/.julia/dev/TrussMorph/test/truss_json/2D_truss.json"
fp1 = "/Users/yijiangh/.julia/dev/TrussMorph/test/truss_json/2D_truss_1.json"
t0 = parse_truss_json(fp0)
t1 = parse_truss_json(fp1)

node_dof = 3
full_node_dof = 3
load_fp = "/Users/yijiangh/.julia/dev/TrussMorph/test/truss_json/2D_truss_load_case.json"
load = parse_load_json(load_fp, node_dof)

morph_path = tm.compute_morph_path(t0, t1, load, node_dof, full_node_dof, path_disc=5)

scene = Scene()
tm.draw_load!(scene, t0, load, load_scale=0.05)

for i=1:length(morph_path)
    tm.draw_truss!(scene, morph_path[i], t0.T, t0.S, r, supp_scale=0.2)
    display(scene)
end
# tm.draw_deformed!(scene, t, U, node_dof)
