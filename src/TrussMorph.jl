
__precompile__()

module TrussMorph

# dependencies
using JSON
using Makie
using LinearAlgebra

# abstract types

# exported APIs
export Truss, MaterialProperties, SectionProperties
export parse_truss_json, parse_load_json, draw_truss!, draw_load!, draw_deformed!
export compute_morph_path

export assemble_global_stiffness_matrix, assemble_load_vector, dof_permutation, get_weight_calculation_fn

include("common_types.jl")
include("json_parser.jl")
include("draw_utils.jl")
include("stiffness.jl")
include("morph.jl")
include("utils.jl")

end # module
