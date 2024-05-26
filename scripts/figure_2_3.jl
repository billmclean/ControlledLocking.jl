import ControlledLocking: optimal_α
import ControlledLocking.Example2: numerical_solution, essential_bcs,
                                   deformed_mesh, selected_displacement,
                                   geometry_filename, diam_Ω, λ, μ, g, u2_A
import SimpleFiniteElements.Elasticity: visualise_soln
import SimpleFiniteElements.Utils: gmsh2pyplot
using SimpleFiniteElements
using PyPlot
using Printf

@printf("\nCooks membrane benchmark problem.\n")
@printf("\tλ = %g, μ = %g\n", λ, μ)
path = joinpath("..", "spatial_domains", geometry_filename)
gmodel = GeometryModel(path)
hmax = 4.0
num_solutions = 5
mesh = FEMesh(gmodel, hmax; order=1, save_msh_file=false,
              refinements=num_solutions-1, verbosity=2)

k = 2
dof = DegreesOfFreedom(mesh[k], essential_bcs)
std_u1h, std_u2h = numerical_solution(λ, μ, 1.0, g, dof)
h = max_elt_diameter(mesh[k])
α = optimal_α(h, λ, diam_Ω)
ctrl_u1h, ctrl_u2h = numerical_solution(λ, μ, α, g, dof)

figure(1, figsize=(11,5))
std_dof = deformed_mesh(std_u1h, std_u2h, dof)
std_x, std_y, std_triangles = gmsh2pyplot(std_dof)
ctrl_dof = deformed_mesh(ctrl_u1h, ctrl_u2h, dof)
ctrl_x, ctrl_y, ctrl_triangles = gmsh2pyplot(ctrl_dof)
x = [0, 48, 48, 0, 0]
y = [0, 44, 60, 44, 0]
subplot(1, 2, 1)
plot(x, y, "--")
triplot(std_x, std_y, std_triangles, linewidth=0.5)
axis("equal")
text(30, 10, "α = 1")
subplot(1, 2, 2)
plot(x, y, "--")
triplot(ctrl_x, ctrl_y, ctrl_triangles, linewidth=0.5)
axis("equal")
s = @sprintf("α = %0.3f", α)
text(30, 10, s)

h = max_elt_diameter(mesh[k])
Nₕ = 2 * dof.num_free
@printf("Showing deformed meshes with h = %g and \
Nₕ = %d degrees of freedom.\n", h, Nₕ)
savefig("deformed_meshes_$Nₕ.pdf")

Nₕ = Vector{Float64}(undef, num_solutions)
std_u1_A = similar(Nₕ)
std_u2_A = similar(Nₕ)
ctrl_u1_A= similar(Nₕ)
ctrl_u2_A = similar(Nₕ)
for k = 1:num_solutions
    global Nₕ, std_u1_A, std_u2_A, ctrl_u1_A, ctrl_u2_A
    local dof, h, std_u1h, std_u2h, α, ctrl_u1h, ctrl_u2h
    dof = DegreesOfFreedom(mesh[k], essential_bcs)
    Nₕ[k] = dof.num_free
    std_u1h, std_u2h = numerical_solution(λ, μ, 1.0, g, dof)
    std_u1_A[k], std_u2_A[k] = selected_displacement(std_u1h, std_u2h, dof)
    h = max_elt_diameter(mesh[k])
    α = optimal_α(h, λ, diam_Ω)
    ctrl_u1h, ctrl_u2h = numerical_solution(λ, μ, α, g, dof)
    ctrl_u1_A[k], ctrl_u2_A[k] = selected_displacement(ctrl_u1h, ctrl_u2h, dof)
end

figure(2)
semilogx(Nₕ, ctrl_u2_A, "o-", 
	 Nₕ, std_u2_A, "x-", 
	[Nₕ[1], Nₕ[end]], [u2_A, u2_A], "k-")
legend((L"$\alpha=\alpha_*(\lambda,h)$", L"$\alpha=1$ "))
grid(true)
xlabel(L"$N_h$")
ylabel(L"$u_2(A)$")
savefig("benchmark.pdf")

