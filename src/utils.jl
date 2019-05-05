using StaticArrays
using LinearAlgebra

function local_coord_rotation(u::SVector{2, Float64}, v::SVector{2, Float64})
    local L = norm(u - v)
    c = (u - v) ./ L
    R::SMatrix{2, 2, Float64} = zeros(2, 2)
    R[1,1] = c[1]
    R[1,2] = c[2]
    R[2,1] = -c[2]
    R[2,2] = c[1]
    return R
end

function compute_round_section_properties(radius::Float64)
    @assert(radius > 0)
    local A = pi * radius^2
    local Jx = 0.5 * pi * radius^4
    local Iy = pi * radius^4 / 4
    sp = SectionProperties(A, Jx, Iy, Iy) # assuming Iz = Iy
    return sp
end
