using TrussMorph
using Test
using Makie

fp = "/Users/yijiangh/.julia/dev/TrussMorph/test/truss_json/2D_truss.json"
truss = parse_truss_json(fp)
load_fp = "/Users/yijiangh/.julia/dev/TrussMorph/test/truss_json/2D_truss_load_case.json"

load = parse_load_json(load_fp)
A = ones(size(truss.T,1))

scene = Scene()
draw_truss!(scene, truss, A, supp_scale=0.2)
draw_load!(scene, truss, load, load_scale=0.05)

# display(scene)

# @test
