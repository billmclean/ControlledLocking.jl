module Example2

import SimpleFiniteElements: DegreesOfFreedom
import SimpleFiniteElements.Elasticity: ∫g_dot_v!, ∫∫λ_div_u_div_v!,
                                       ∫∫2μ_εu_εv!, elasticity_soln, error_norms
import SimpleFiniteElements.FEM: prepare_assembly
import SimpleFiniteElements.Utils: barycentric!
import SimpleFiniteElements.MeshGen: elt_node_coord!
import StaticArrays: SA

const diam_Ω = hypot(48, 60)
const geometry_filename = "cooks_membrane.geo"
gD = SA[0.0, 0.0]
const essential_bcs = [("Left", gD)]

const E = 1.12499998125
const ν = 0.499999975
const λ = E * ν / ( (1 + ν) * (1 - 2ν) )
const μ = E / ( 2 * (1 + ν) )
const g = 1/16
const u2_A = 16.442

function numerical_solution(λ::Float64, μ::Float64, α::Float64, c::Float64,
	dof::DegreesOfFreedom)
    bilinear_forms = Dict("Omega" => [(∫∫λ_div_u_div_v!, λ^α),
                                      (∫∫2μ_εu_εv!, μ)])
    gN = SA[0.0, c]
    linear_funcs = Dict("Right" => (∫g_dot_v!, gN))
    u1h, u2h = elasticity_soln(dof, bilinear_forms, linear_funcs)
    return u1h, u2h
end

function divergence(u1h::Vector{Float64}, u2h::Vector{Float64}, 
                    dof::DegreesOfFreedom, name::String)
    mesh = dof.mesh
    ψ = Vector{Float64}(undef, 3)
    b = Matrix{Float64}(undef, 2, 3)
    cntrd = Vector{Float64}(undef, 2)
    elt_type = mesh.elt_type_in[name]
    if elt_type ≠ 2
	error("Elements in $name must be of type Triangle 3")
    end
    each_elt_index, _, _, elt_global_dof, _,
       coord = prepare_assembly(name, dof, 1)
    num_elts = length(each_elt_index)
    div_uh = Vector{Float64}(undef, num_elts)
    k = 1
    L2norm = 0.0
    total_area = 0.0
    for l in each_elt_index
        elt_node_coord!(coord, elt_global_dof, name, l, dof)
        area = barycentric!(b, cntrd, coord)
	total_area += area
	∂₁u₁ = 0.0
	∂₂u₂ = 0.0
        for j = eachindex(elt_global_dof)
            r = elt_global_dof[j]
	    ∂₁u₁ += u1h[r] * b[1,j]
	    ∂₂u₂ += u2h[r] * b[2,j]
        end
	div_uh[k] = ∂₁u₁ + ∂₂u₂
	L2norm += area * div_uh[k]^2
	k += 1
    end
    L2norm = sqrt(L2norm/total_area)
    return div_uh, L2norm, total_area
end

function selected_displacement(u1h::Vector{Float64}, u2h::Vector{Float64},
	                            dof::DegreesOfFreedom)
    A = [48, 52]
    mesh = dof.mesh
    elt_node_tags = mesh.elt_node_tags_in["Right"]
    nelts = size(elt_node_tags, 2)
    if iseven(nelts)
	k = nelts ÷ 2
	nd = elt_node_tags[2,k]
	P = mesh.coord[:,nd]
	if maximum(abs, A - P) > eps(Float64)
	    error("Failed to find the correct node")
	end
	k = dof.node_tag_to_dof[nd]
	u1_A = u1h[k]
	u2_A = u2h[k]
    else
    end
    return u1_A, u2_A
end

function deformed_mesh(u1h::Vector{Float64}, u2h::Vector{Float64}, 
	               dof::DegreesOfFreedom)
    new_dof = deepcopy(dof)
    new_coord = new_dof.mesh.coord
    for k = 1:dof.num_free
	nd = dof.node_tag[k]
	new_coord[1,nd] += u1h[k]
	new_coord[2,nd] += u2h[k]
    end
    return new_dof
end

end # module
