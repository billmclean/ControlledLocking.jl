import ControlledLocking: optimal_α
import ControlledLocking.Example2: numerical_solution, essential_bcs,
                                   Lame_parameters, divergence,
				   geometry_filename, diam_Ω, E, ν, λ, μ, g
using SimpleFiniteElements
using Printf

@printf("Cooks membrane problem with\n")
@printf("\tE = %g\n", E)
@printf("\tν = %g\n", ν)
@printf("\tλ = %g\n", λ)
@printf("\tμ = %g\n", μ)
@printf("Traction constant g = %g\n", g)

path = joinpath("..", "spatial_domains", geometry_filename)
gmodel = GeometryModel(path)
hmax = 4.0
nrows = 6
mesh = FEMesh(gmodel, hmax; order=1, save_msh_file=false,
              refinements=nrows-1, verbosity=2)
@printf("\n%6s  %5s  %5s  %8s    %5s  %8s\n\n", 
	"Nₕ", "h", "α", "||∇⋅uₕ||", "α", "||∇⋅uₕ||")
for row = 1:nrows
    local Nₕ
    dof = DegreesOfFreedom(mesh[row], essential_bcs)
    std_u1h, std_u2h = numerical_solution(λ, μ, 1.0, g, dof)
    std_div_uh, std_L2norm, area = divergence(std_u1h, std_u2h, dof, "Omega")
    h = max_elt_diameter(mesh[row])
    α = optimal_α(h, λ, diam_Ω)
    ctrl_u1h, ctrl_u2h = numerical_solution(λ, μ, α, g, dof)
    ctrl_div_uh, ctrl_L2norm, _ = divergence(ctrl_u1h, ctrl_u2h, dof, "Omega")
    Nₕ = 2 * dof.num_free
    @printf("%6d  %5.3f  %5.3f  %8.2e    %5.3f  %8.2e\n", 
	    Nₕ, h, 1.0, std_L2norm, α, ctrl_L2norm)
end

