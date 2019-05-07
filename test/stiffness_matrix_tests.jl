using TrussMorph
tm = TrussMorph
using LinearAlgebra: cholesky, I
using Test
using Makie

fp = "/Users/yijiangh/.julia/dev/TrussMorph/test/truss_json/2D_truss_1.json"
t = tm.parse_truss_json(fp)
load_fp = "/Users/yijiangh/.julia/dev/TrussMorph/test/truss_json/2D_truss_load_case.json"

# truss
node_dof = 2
full_node_dof = 3

load = tm.parse_load_json(load_fp, node_dof)

n_v = size(t.X,1)
n_e = size(t.T,1)
# initial cross secs
A = ones(size(t.T,1))* pi * 0.01^2
sys_dof = node_dof * n_v

box_range = 6
t.X[2:3,:] += hcat([-1, -1] .+ 2 * rand(Float64,2), [1, 1]*(-box_range/2) .+ box_range * rand(Float64,2))

K, KR_es, id_map = tm.assemble_global_stiffness_matrix(t.X, t.T, A, t.mp.E, node_dof, full_node_dof)
F = tm.assemble_load_vector(load, n_v, node_dof)

perm_vec, perm, n_dof_free = tm.dof_permutation(t.S, n_v, node_dof)
K_perm = perm * K * perm'
K_mm = K_perm[1:n_dof_free, 1:n_dof_free]

F_perm = perm * F
F_m = Array(F_perm[1:n_dof_free])

# K_mm = (K_mm + K_mm')/2
Kpf = cholesky(K_mm)
U_m = Kpf\F_m
# U_I, U_V = findnz(U_m)
U_perm = vcat(U_m, zeros(sys_dof - n_dof_free))
U = perm' * U_perm

weight = tm.weight_calculation(t.X, A, t.T, F_m, perm, n_dof_free, node_dof, full_node_dof, t.mp.E)

scene = Scene()
tm.draw_truss!(scene, t, A, supp_scale=0.2)
tm.draw_load!(scene, t, load, load_scale=0.05)
# tm.draw_deformed!(scene, t, U, node_dof)

# @test output_string(:x) == "x"
# @test_throws MethodError output_string("MySymbol")
