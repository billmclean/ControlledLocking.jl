module Example1

import SimpleFiniteElements: DegreesOfFreedom
import SimpleFiniteElements.Elasticity: ∫∫f_dot_v!, ∫∫λ_div_u_div_v!,
                                       ∫∫2μ_εu_εv!, elasticity_soln, error_norms
import StaticArrays: SA

const diam_Ω = sqrt(2)
const geometry_filename = "square.geo"

function exact_u(x, y, λ)
    common_term = sin(x) * sin(y) / λ
    u1 = ( cos(2x) - 1 ) * sin(2y) + common_term
    u2 = ( 1 - cos(2y) ) * sin(2x) + common_term
    return SA[u1, u2]
end

function ∇u(x, y, λ)
    sx, cx = sincos(x)
    sy, cy = sincos(y)
    s2x, c2x = sincos(2x)
    s2y, c2y = sincos(2y)
    s2x_s2y = s2x * s2y
    cx_sy = cx * sy
    sx_cy = sx * cy
    ∂₁u₁ = -2s2x_s2y        + cx_sy / λ
    ∂₂u₁ = 2(c2x - 1) * c2y + sx_cy / λ
    ∂₁u₂ = 2(1 - c2y) * c2x + cx_sy / λ
    ∂₂u₂ =  2s2x_s2y        + sx_cy / λ
    return SA[ ∂₁u₁  ∂₁u₂
               ∂₂u₁  ∂₂u₂ ]
end

function ϕ(x, y)
    return 4sin(2y) * (2cos(2x) - 1)
end

function f(x, y, λ)
    ψ1 = 2sin(x) * sin(y) - cos(x+y) 
    ψ3 = cos(x+y)
    f1 =  ϕ(x, y) - ψ3 + ψ1 / λ
    f2 = -ϕ(y, x) - ψ3 + ψ1 / λ
    return SA[f1, f2]
end

gD = SA[0.0, 0.0]
const essential_bcs = [("Top", gD), ("Bottom", gD), ("Left", gD), ("Right", gD)]
const μ = 1.0

function numerical_solution(λ::Float64, λₕ::Float64, dof::DegreesOfFreedom;
                            μ=1.0)
    bilinear_forms = Dict("Omega" => [(∫∫λ_div_u_div_v!, λₕ),
                                      (∫∫2μ_εu_εv!, μ)])
    linear_funcs = Dict("Omega" => (∫∫f_dot_v!, f, λ))
    u1h, u2h = elasticity_soln(dof, bilinear_forms, linear_funcs)
    return u1h, u2h, dof
end

function numerical_errors(u1h, u2h, dof::DegreesOfFreedom, λ::Float64;
                          quadrature_level=2)
    L2err, H1err = error_norms(u1h, u2h, exact_u, ∇u, dof, 
                               quadrature_level, λ)
    return L2err, H1err
end

function numerical_errors(λ::Float64, λₕ::Float64, dof::DegreesOfFreedom;
                          quadrature_level=2)
    u1h, u2h, dof = numerical_solution(λ, λₕ, dof)
    L2err, H1err = numerical_errors(u1h, u2h, dof, λ; 
				    quadrature_level=quadrature_level)
    return L2err, H1err
end

end # module
