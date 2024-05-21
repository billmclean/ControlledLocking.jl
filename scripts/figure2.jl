import ControlledLocking: optimal_α
import ControlledLocking.Example2: numerical_solution, essential_bcs, 
                                   Lame_parameters, divergence,
				   geometry_filename, diam_Ω, λ, μ, g
import SimpleFiniteElements.Elasticity: visualise_soln
import SimpleFiniteElements.Utils: gmsh2pyplot
using SimpleFiniteElements
using PyPlot
using Printf

@printf("Cooks membrane problem with\n")
@printf("\tλ = %g\n", λ)
@printf("\tμ = %g\n", μ)
@printf("Traction constant g = %g\n", g)

path = joinpath("..", "spatial_domains", geometry_filename)
gmodel = GeometryModel(path)
hmax = 2.0
mesh = FEMesh(gmodel, hmax; order=1, save_msh_file=true,
              refinements=0, verbosity=2)
dof = DegreesOfFreedom(mesh, essential_bcs)
@printf("Number of degrees of freedom = %d\n", 2dof.num_free)

std_u1h, std_u2h = numerical_solution(λ, μ, 1.0, c, dof)
std_div_uh, _, _ = divergence(std_u1h, std_u2h, dof, "Omega")
h = max_elt_diameter(mesh)
α = optimal_α(h, λ, diam_Ω)
@printf("α = %g\n", α)
ctrl_u1h, ctrl_u2h = numerical_solution(λ, μ, α, c, dof)
ctrl_div_uh, _, _ = divergence(ctrl_u1h, ctrl_u2h, dof, "Omega")

x, y, triangles = gmsh2pyplot(dof)
figure(1, figsize=(11,5))
subplot(1, 2, 1)
tripcolor(x, y, triangles, std_div_uh)
text(30, 10, "α = 1")
colorbar()
axis("equal")
subplot(1, 2, 2)
tripcolor(x, y, triangles, ctrl_div_uh)
s = @sprintf("α = %g", α)
text(30, 10, s)
colorbar()
axis("equal")
filename = @sprintf("divergence_%d.pdf", 2dof.num_free)
savefig(filename)

scale = 5.0
visualise_soln(dof,  std_u1h,  std_u2h, scale, 3)
visualise_soln(dof, ctrl_u1h, ctrl_u2h, scale, 5)
