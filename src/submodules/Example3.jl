module Example3

import SimpleFiniteElements: DegreesOfFreedom
import SimpleFiniteElements.Elasticity: ∫∫f_dot_v!, ∫∫λ_div_u_div_v!,
                                       ∫∫2μ_εu_εv!, elasticity_soln, error_norms
import StaticArrays: SA

const diam_Ω = sqrt(2) * π
const geometry_filename = "square.geo"

α = 0.5
β = 1.0
τ(x, y, α) = 1 + α * sin(2x)
λ(x, y, α, Λ) = Λ * τ(x, y, α)
μ(x, y, β) = 1 + β * (x + y)
μ_bar = 1 + β

function exact_u(x, y, Λ)
    common_term = sin(x) * sin(y) / Λ
    u1 = ( cos(2x) - 1 ) * sin(2y) + common_term
    u2 = ( 1 - cos(2y) ) * sin(2x) + common_term
    return SA[u1, u2]
end

function ∇u(x, y, Λ)
    sx, cx = sincos(x)
    sy, cy = sincos(y)
    s2x, c2x = sincos(2x)
    s2y, c2y = sincos(2y)
    s2x_s2y = s2x * s2y
    cx_sy = cx * sy
    sx_cy = sx * cy
    ∂₁u₁ = -2s2x_s2y        + cx_sy / Λ
    ∂₂u₁ = 2(c2x - 1) * c2y + sx_cy / Λ
    ∂₁u₂ = 2(1 - c2y) * c2x + cx_sy / Λ
    ∂₂u₂ =  2s2x_s2y        + sx_cy / Λ
    return SA[ ∂₁u₁  ∂₁u₂
	       ∂₂u₁  ∂₂u₂ ]
end


function ϕ(x, y, β)
    sx, cx = sincos(2x)
    sy, cy = sincos(2y)
    return 2β * ( 2sx * sy - cx + cy ) + 4μ(x, y, β) * sy * (2cx - 1)
end

function f(x, y, α, β, Λ)
    ψ1 = μ(x, y, β) * ( 2sin(x) * sin(y) - cos(x+y) ) - β * sin(x+y)
    ψ2(x, y) = 2β * cos(x) * sin(y)
    ψ3 = τ(x, y, α) * cos(x+y)
    f1 =  ϕ(x, y, β) - ψ3 - 2α * cos(2x) * sin(x+y) + (ψ1 -  ψ2(x, y)) / Λ
    f2 = -ϕ(y, x, β) - ψ3 + (ψ1 - ψ2(y, x)) / Λ
    return SA[f1, f2]
end

gD = SA[0.0, 0.0]
const essential_bcs = [("Top", gD), ("Bottom", gD), ("Left", gD), ("Right", gD)]

function numerical_solution(Λ::Float64, Λₕ::Float64, dof::DegreesOfFreedom, μ::Function)
    #bilinear_forms = Dict("Omega" => [(∫∫λ_div_u_div_v!, λₕ),
    #                                  (∫∫2μ_εu_εv!, μ)])
    #linear_funcs = Dict("Omega" => (∫∫f_dot_v!, f, λ))
    bilinear_forms = Dict("Omega" => [(∫∫λ_div_u_div_v!, λ, α, Λₕ),
                                      (∫∫2μ_εu_εv!, μ, β)])
    linear_funcs = Dict("Omega" => (∫∫f_dot_v!, f, α, β, Λ))
    u1h, u2h = elasticity_soln(dof, bilinear_forms, linear_funcs)
    return u1h, u2h, dof
end

function numerical_errors(u1h, u2h, dof::DegreesOfFreedom, Λ::Float64;
                          quadrature_level=2)
    L2err, H1err = error_norms(u1h, u2h, exact_u, ∇u, dof, 
                               quadrature_level, Λ)
    return L2err, H1err
end

function numerical_errors(Λ::Float64, Λₕ::Float64, dof::DegreesOfFreedom;
                          quadrature_level=2)
    u1h, u2h, dof = numerical_solution(Λ, Λₕ, dof, μ)
    L2err, H1err = numerical_errors(u1h, u2h, dof, Λ; 
				    quadrature_level=quadrature_level)
    return L2err, H1err
end

end # module
