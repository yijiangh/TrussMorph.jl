using TrussMorph
using Test

fp = "/Users/yijiangh/.julia/dev/TrussMorph/test/truss_json/2D_truss.json"
truss = parse_truss_json(fp)

load_fp = "/Users/yijiangh/.julia/dev/TrussMorph/test/truss_json/2D_truss_load_case.json"
load = parse_load_json(load_fp)

@test size(truss.X,2) == 2
@test size(truss.T,2) == 2
@test size(truss.S,2) == 4
@test size(load, 2) == 4
# @test_throws MethodError parse_truss_json(fp)
