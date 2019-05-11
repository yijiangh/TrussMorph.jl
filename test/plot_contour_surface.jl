using TrussMorph
tm = TrussMorph
using Test
using Makie
using Colors

node_dof = 3 # 2 for truss model, 3 for frames
full_node_dof = 3 # fixed for 2D cases

dir = "/Users/yijiangh/Dropbox (MIT)/Course_Work/6.838_2019_spring/final_project/code/gh_validation"
result_file_dir = joinpath(pwd(),"test","results")

st_file_name  = "2D_truss_1.json"
end_file_name  = "2D_truss.json"
fp0 = joinpath(dir, st_file_name)
fp1 = joinpath(dir, end_file_name)
load_fp = joinpath(dir, "2D_truss_load_case.json")
design_var_ids = [4, 6]

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

var_chuck = length(design_var_ids)::Int
X0_var = reshape(t0.X', prod(size(t0.X)))[design_var_ids]
X1_var = reshape(t1.X', prod(size(t1.X)))[design_var_ids]

weight_fn = tm.get_weight_calculation_fn(t0.X, r, t0.T, F_m, perm, n_dof_free,
node_dof, full_node_dof, t0.mp, design_var_ids)

y1 = -2:0.01:2h
y2 = -2:0.01:2

f(x, y) = log(weight_fn([x,y]))

sc = Scene()
p1 = surface!(sc, y1, y2, f)
# p2 = contour3d!(sc, y1, y2, (x, y) -> f(x,y), levels = 15, linewidth = 3)

display(sc)
