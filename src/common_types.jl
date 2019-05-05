
mutable struct MaterialProperties
    E::Float64 # Young's modulus
    G::Float64 # Shear Modulus
    μ::Float64 # Poisson ratio
    ρ::Float64 # density
end

mutable struct SectionProperties
    A::Float64
    Jx::Float64
    Iy::Float64
    Iz::Float64
end

mutable struct Truss
    X::Matrix{Float64} # nodes x 2 matrix
    T::Matrix{Int}     # elements x 2 matrix
    S::Matrix{Int}     # fixities x 4 matrix
    mp::MaterialProperties
end
