module ControlledLocking

include("submodules/Example1.jl")
include("submodules/Example2.jl")

function optimal_α(h::Float64, λ::Float64, diam_Ω::Float64)
    α = min(1, log(diam_Ω/h) / log(λ))
end

function reduced_λ(h::Float64, λ::Float64, diam_Ω::Float64)
    λₕ = λ / ( 1 + λ * h / diam_Ω )
end

end # module ControlledLocking
