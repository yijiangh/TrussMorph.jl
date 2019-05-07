
__precompile__()

module TrussMorph

# dependencies
using JSON
using Makie
using LinearAlgebra

# abstract types

# exported APIs
export parse_truss_json, parse_load_json, draw_truss!, draw_load!, draw_deformed!

# export assemble_global_stiffness_matrix

include("common_types.jl")
include("json_parser.jl")
include("draw_utils.jl")

end # module
