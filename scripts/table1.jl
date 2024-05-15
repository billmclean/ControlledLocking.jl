import ControlledLocking.Example1: numerical_errors, essential_bcs
using SimpleFiniteElements
using Printf

path = joinpath("..", "spatial_domains", "square.geo")
gmodel = GeometryModel(path)
hmax = 0.2
nrows = 6
mesh = FEMesh(gmodel, hmax, order=1, save_msh_file=false,
              refinements=nrows-1, verbosity=2)

std_L2err = Vector{Float64}(undef, nrows)
std_H1err = similar(std_L2err)
ctrl_L2err = similar(std_L2err)
ctrl_H1err = similar(std_L2err)
Nₕ = Vector{Int64}(undef, nrows)

λ = 10000.0
@printf("Conforming elements, λ = %g\n", λ)

@printf("\nL2 errors:\n")
@printf("%6s  %6s  %5s  %8s  %5s  %10s  %5s  %8s\n\n", 
        "h", "DoF", "α", "Standard", "rate", "Controlled", "rate", "secs")

for row in eachindex(std_L2err)
    start = time()
    dof = DegreesOfFreedom(mesh[row], essential_bcs)
    std_L2err[row], std_H1err[row] = numerical_errors(λ, 1.0, dof)
    h = max_elt_diameter(mesh[row])
    α = min(1, log(1/h) / log(λ))
    Nₕ[row] = dof.num_free
    ctrl_L2err[row], ctrl_H1err[row] = numerical_errors(λ, α, dof)
    elapsed = time() - start
    if row == 1
        @printf("%6.3f  %6d  %5.3f  %8.2e  %5s  %10.2e  %5s  %8.3f\n", 
                h, Nₕ[row], α, std_L2err[row], "", ctrl_L2err[row], "", elapsed)
    else
        std_rate = log2(std_L2err[row-1] / std_L2err[row])
        ctrl_rate = log2(ctrl_L2err[row-1] / ctrl_L2err[row])
        @printf("%6.3f  %6d  %5.3f  %8.2e  %5.3f  %10.2e  %5.3f  %8.3f\n", 
                h, Nₕ[row], α, std_L2err[row], std_rate, 
                ctrl_L2err[row], ctrl_rate, elapsed)
    end
end

@printf("\nH1 errors:\n")
@printf("%6s  %6s  %8s  %5s  %10s  %5s\n\n", 
        "h", "DoF", "Standard", "rate", "Controlled", "rate")

for row in eachindex(std_L2err)
    start = time()
    h = max_elt_diameter(mesh[row])
    α = min(1, log(1/h) / log(λ))
    if row == 1
        @printf("%6.3f  %6d  %8.2e  %5s  %10.2e  %5s\n", h, Nₕ[row], 
                std_H1err[row], "", ctrl_H1err[row], "")
    else
        std_rate = log2(std_H1err[row-1] / std_H1err[row])
        ctrl_rate = log2(ctrl_H1err[row-1] / ctrl_H1err[row])
        @printf("%6.3f  %6d  %8.2e  %5.3f  %10.2e  %5.3f\n",
                h, Nₕ[row], std_H1err[row], std_rate, 
                ctrl_H1err[row], ctrl_rate)
    end
end
