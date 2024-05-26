import ControlledLocking: optimal_α
import ControlledLocking.Example1: numerical_errors, essential_bcs,
				   geometry_filename, diam_Ω
using SimpleFiniteElements
using Printf

path = joinpath("..", "spatial_domains", geometry_filename)
gmodel = GeometryModel(path)
hmax = 0.4
nrows = 6
mesh = FEMesh(gmodel, hmax, order=1, save_msh_file=false,
              refinements=nrows-1, verbosity=2)
dof = Vector{DegreesOfFreedom}(undef, nrows)
h = Vector{Float64}(undef, nrows)
for row in eachindex(mesh)
    dof[row] = DegreesOfFreedom(mesh[row], essential_bcs)
    h[row] = max_elt_diameter(mesh[row])
end

λ = [1_000.0, 10_000.0, 100_000.0]
nblocks = length(λ)
block_method = [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2)]

L2err = Array{Float64}(undef, nrows, nblocks, 2)
H1err = similar(L2err)
α = Array{Float64}(undef, nrows, nblocks)

ncases = 2 * nblocks
@printf("Using %d threads to handle %d cases.\n", 
	Threads.nthreads(), ncases)
start = time()
@Threads.threads for k = 1:ncases
    block, method = block_method[k]
    for row = 1:nrows
	if method == 1
	    L2err[row,block,1], 
	    H1err[row,block,1] = numerical_errors(
                                  λ[block], 1.0, dof[row]) 
	else
	    α[row,block] = optimal_α(h[row], λ[block], diam_Ω)
	    L2err[row,block,2], 
	    H1err[row,block,2] = numerical_errors(
                                  λ[block], α[row,block], dof[row]) 
	end
    end
end
elapsed = time() - start
@printf("PDE solves took %g seconds\n", elapsed)

function print_table(io::IO, err::Array{Float64}, 
	             block::Int64, title::String)
    @printf(io, "\n%70s\n", repeat('-', 70))
    @printf(io, "\n%s\n\n", title)
    @printf(io, "%12s      %8s           %10s\n", "", "Standard", "Controlled")
    @printf(io, "%6s  %5s  %8s  %5s    %8s  %5s  %5s\n\n",
	    "Nₕ", "h", "error", "rate", "error", "rate", "α")
    for row = 1:nrows
        Nₕ = 2 * dof[row].num_free
	@printf(io, "%6d& %5.3f& ", Nₕ, h[row])
	if row == 1
	    @printf(io, "%8.2e& %5s&   %8.2e& %5s& %5.3f ", 
	            err[row,block,1], "", err[row,block,2], "", α[row,block])
	else
	    rate1 = log2(err[row-1,block,1] / err[row,block,1])
	    rate2 = log2(err[row-1,block,2] / err[row,block,2])
            @printf(io, "%8.2e& %5.3f&   %8.2e& %5.3f& %5.3f ", 
                    err[row,block,1], rate1, 
		    err[row,block,2], rate2, α[row,block])
	end
	@printf(io, "\\\\\n")
    end
end

output_file = "table1_output.txt"
@printf("Writing output to %s\n", output_file)
open(output_file, "w") do io
    @printf(io, "Output from table1.jl\n")
    for block = 1:nblocks
        title = @sprintf("L2 errors, λ = %g", λ[block])
        print_table(io, L2err, block, title)
    end
    for block = 1:nblocks
        title = @sprintf("H1 errors, λ = %g", λ[block])
        print_table(io, H1err, block, title)
    end
end

