using TrussMorph
tm = TrussMorph
using Test
using Makie

fp = "/Users/yijiangh/.julia/dev/TrussMorph/test/truss_json/2D_truss_1.json"
t = tm.parse_truss_json(fp)
load_fp = "/Users/yijiangh/.julia/dev/TrussMorph/test/truss_json/2D_truss_load_case.json"

# truss
node_dof = 3 # 2 for truss, 3 for frame
full_node_dof = 3

load = tm.parse_load_json(load_fp, node_dof)

n_v = size(t.X,1)
n_e = size(t.T,1)

# TODO: perturb this and see resulting weight variance
# initial cross secs radius
r = ones(size(t.T,1)) * 0.01

# pertubing the truss randomly within a box
iter = 10
# scene = Scene()

st_t = time()
for i=1:iter
    box_range = 6
    X = copy(t.X)
    X[2:3,:] += hcat([-1, -1] .+ 2 * rand(Float64,2), [1, 1]*(-box_range/2) .+ box_range * rand(Float64,2))

    F = tm.assemble_load_vector(load, n_v, node_dof)
    perm_vec, perm, n_dof_free = tm.dof_permutation(t.S, n_v, node_dof)

    F_perm = perm * F
    F_m = Array(F_perm[1:n_dof_free])

    weight = tm.weight_calculation(X, r, t.T, F_m, perm, n_dof_free, node_dof, full_node_dof, t.mp)
    @show weight
end
println("avg time: ", (time() - st_t) / iter)

    # tm.draw_truss!(scene, t, A, supp_scale=0.2)
    # tm.draw_load!(scene, t, load, load_scale=0.05)
    # tm.draw_deformed!(scene, t, U, node_dof)
