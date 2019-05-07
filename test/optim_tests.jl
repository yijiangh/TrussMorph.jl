using Optim
using LinearAlgebra: I, det, tr

# Rosenbrock
f(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
x0 = [0.0, 0.0]
optimize(f, x0, LBFGS(), autodiff = :forward)

# analytical gradient
function g!(G, x)
    G[1] = -2.0 * (1.0 - x[1]) - 400.0 * (x[2] - x[1]^2) * x[1]
    G[2] = 200.0 * (x[2] - x[1]^2)
end
optimize(f, g!, x0)

# analytical Hessian
function h!(H, x)
    H[1, 1] = 2.0 - 400.0 * x[2] + 1200.0 * x[1]^2
    H[1, 2] = -400.0 * x[1]
    H[2, 1] = -400.0 * x[1]
    H[2, 2] = 200.0
end
optimize(f, g!, h!, x0, LBFGS())

lower = [1.25, -2.1]
upper = [Inf, Inf]
initial_x = [2.0, 2.0]
# requires using LineSearches
inner_optimizer = GradientDescent()
res = optimize(f, g!, lower, upper, initial_x, Fminbox(inner_optimizer))

summary(results)
Optim.minimizer(results)
Optim.minimum(results)


function fMat(Xm::Matrix{Float64})
    # Xm = reshape(X,(3,3))
    return det((Xm-I) * (Xm-I))
end
X0 = reshape(collect(1:1.0:9.0), (3,3))
optimize(fMat, X0, LBFGS())

function sq(x)
    return x^2
end

function get_fn(a, b)::Function
    function calc_fn(x::Vector{Float64})
        return sum((x .* sq(a) .+ b).^2)
    end
    return calc_fn
end

f = get_fn(2.0,3.0)
optimize(f, [1.0, 2.0], LBFGS())
