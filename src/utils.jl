function compute_round_section_properties(radius::Float64)
    @assert(radius > 0)
    local A = pi * radius^2
    local Jx = 0.5 * pi * radius^4
    local Iy = pi * radius^4 / 4
    sp = SectionProperties(A, Jx, Iy, Iy) # assuming Iz = Iy
    return sp
end
