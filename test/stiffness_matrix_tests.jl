using TrussMorph
using Test

fp = "/Users/yijiangh/.julia/dev/TrussMorph/test/truss_json/2D_truss.json"
t = parse_truss_json(fp)
load_fp = "/Users/yijiangh/.julia/dev/TrussMorph/test/truss_json/2D_truss_load_case.json"

# truss
node_dof = 2
full_node_dof = 3

load = parse_load_json(load_fp, node_dof)

n_v = size(t.X,1)
n_e = size(t.T,1)
# initial cross secs
A = ones(size(t.T,1))
sys_dof = node_dof * n_v

K = assemble_global_stiffness_matrix(t.X, t.T, A, t.mp.E, node_dof, full_node_dof)
F = assemble_load_vector(load, n_v, node_dof)

_, perm, n_dof_free = dof_permutation(t.S, n_v, node_dof)
K_perm = perm * K * perm'
K_mm = K_perm[1:n_dof_free, 1:n_dof_free]

F_perm = perm * F
F_m = Array(F_perm[1:n_dof_free])

U_m = K_mm\F_m
# U_I, U_V = findnz(U_m)
U_perm = vcat(U_m, zeros(sys_dof - n_dof_free))
U = perm' * U_perm

# @test output_string(:MySymbol) == "MySymbol"
# @test output_string(:x) == "x"
# @test_throws MethodError output_string("MySymbol")
