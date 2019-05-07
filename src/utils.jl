# using TrussMorph::SectionProperties

function compute_round_section_properties(radius::Float64)
    # @assert(radius > 0)
    A = pi * radius^2
    Jx = 0.5 * pi * radius^4
    Iy = pi * radius^4 / 4
    return SectionProperties(A, Jx, Iy, Iy) # assuming Iz = Iy
end
