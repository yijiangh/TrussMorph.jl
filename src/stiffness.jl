using TrussMorph

"""
local_stiffness_matrix: compute elemental FRAME stiffness matrix
in the local coordinates.

Return: 6 x 6 matrix K_le (local elemental)

Note: only support 2d right now.
"""
function local_stiffness_matrix(L::Float64, sp::SectionProperties, mp::MaterialProperties)
    E = mp.E
    A = sp.A
    # local G = mp.G
    # local μ = mp.μ

    K_le::Matrix{3,2,Float64} = zeros(2,2)
    K_le[1,1] = 1.0
    K_le[1,2] = -1.0
    K_le[2,1] = 1.0
    K_le[2,2] = -1.0
    K_le .*= E * A / L
    return K_le
end

function assembly_global_stiffness_matrix()


end
