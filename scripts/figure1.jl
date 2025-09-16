import ControlledLocking: reduced_λ
import ControlledLocking.Example1: numerical_errors, essential_bcs,
                                   geometry_filename, diam_Ω
using PyPlot
using SimpleFiniteElements
using Printf

path = joinpath("..", "spatial_domains", geometry_filename)
gmodel = GeometryModel(path)
hmax = 0.4
nrows = 5
mesh = FEMesh(gmodel, hmax, order=1, save_msh_file=false,
              refinements=nrows-1, verbosity=2)
dof = DegreesOfFreedom(mesh[nrows], essential_bcs)
h = max_elt_diameter(mesh[nrows])

λ = 100_000.0
λₕ = range(40, 160, 16)
L = diam_Ω
λₕ_star = λ / sqrt(1 + (λ*h/L)^2)
λₕ_tilde = min(λ, L/h)
λₕ_sharp = reduced_λ(h, λ, L)
Nₕ = 2 * dof.num_free
@printf("λ = %g\n", λ)
@printf("Number of unknowns Nₕ = %d\n", Nₕ)
@printf("Mesh size h = %0.3f\n", h)
@printf("λₕ star  = %0.2f\nλₕ tilde = %0.2f\nλₕ sharp = %0.2f\n",
        λₕ_star, λₕ_tilde, λₕ_sharp)
L2err = zeros(length(λₕ))
H1err = zeros(length(λₕ))
for k in eachindex(λₕ)
    L2err[k], H1err[k] = numerical_errors(λ, λₕ[k], dof)
end

figure(1)
plot(λₕ, L2err, λₕ, H1err)
xmin, xmax, ymin, ymax = axis()
plot([λₕ_sharp, λₕ_sharp], [ymin, ymax], ":")
legend((L"$L^2$-error", L"$H^1$-error"))
grid(true)
xlabel(L"$\lambda_h$")
savefig("error_vs_lambda_h.pdf")

long_λₕ = logrange(10, 10_000, 16)
long_L2err = zeros(length(λₕ))
long_H1err = zeros(length(λₕ))
for k in eachindex(λₕ)
    long_L2err[k], long_H1err[k] = numerical_errors(λ, long_λₕ[k], dof)
end

figure(2)
semilogx(long_λₕ, long_L2err, long_λₕ, long_H1err)
xmin, xmax, ymin, ymax = axis()
plot([λₕ_sharp, λₕ_sharp], [ymin, ymax], ":")
legend((L"$L^2$-error", L"$H^1$-error"))
grid(true)
xlabel(L"$\lambda_h$")
savefig("long_error_vs_lambda_h.pdf")
