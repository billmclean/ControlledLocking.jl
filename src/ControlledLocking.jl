module ControlledLocking

include("submodules/Example1.jl")
include("submodules/Example2.jl")
include("submodules/Example3.jl")

function optimal_α(h::Float64, λ::Float64, diam_Ω::Float64)
    α = min(1, log(diam_Ω/h) / log(λ))
end

function reduced_λ(h::Float64, λ::Float64, diam_Ω::Float64)
    λₕ = λ / ( 1 + λ * h / diam_Ω )
end

function reduced_Λ(h::Float64, Λ::Float64, diam_Ω::Float64, μ_bar::Float64)
    Λ_h = Λ * μ_bar / ( μ_bar + Λ * h / diam_Ω )
end

end # module ControlledLocking
