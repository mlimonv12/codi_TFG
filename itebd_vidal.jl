# iTEBD evolution using vidal notation for the unit cell
using ITensors, ITensorMPS
using Plots
gr()

#include("functions.jl")
anim = Animation()
anim2 = Animation()



# Initializes an MPS as a Vector{ITensor}. Allows for iMPS selection
function init_iMPS2(N::Int64)
    
    sites = siteinds("Electron", N)
    links = [Index(bdim, "link-$i") for i in 0:2*N]
    psi = [ITensor() for _ in 1:(2*N + 1)]

    for i in 1:2*N
        if isodd(i) # Γ
            psi[i] = randomITensor(links[i], sites[i], )
        else # Λ
            psi[i] = randomITensor(links[i], links[i+1])
        end
    end

    psi[N] = psi[1] # Facillitates operations later on

    return psi, sites, links
end




# Generates an iMPS ITensor in Vidal form with a given initial state. 1:empty, 2:up, 3:dn, 4:updn
function init_iMPS(state::Vector; bdim=16)
    
    references = Dict(0=>"Emp", 1=>"Up", 2=>"Dn", 3=>"UpDn")
    
    # Initializes iMPS with two-site unit cell 
    sites = siteinds("Electron", 2)
    links = [Index(bdim, "link-$i") for i in 0:4]
    psi = [ITensor() for _ in 1:4]

    for i in 1:4
        if iseven(i) # Γ
            #T = ITensor(links[i], sites[Int(i/2)], links[i+1])
            #T[links[i]=>bdim, sites[Int(i/2)]=>references[state[Int(i/2)]], links[i+1]=>bdim] = 1.0
            #psi[i] = T
            psi[i] = randomITensor(links[i], sites[Int(i/2)], links[i+1])
        else # Λ
            psi[i] = randomITensor(links[i], links[i+1])
        end
    end

    #psi[5] = psi[1] # Facillitates operations later on

    return psi, sites, links
end


function normalise_unitcell(unitcell::Vector{ITensor})

    unitcell[1] /= norm(unitcell[1])
    unitcell[3] /= norm(unitcell[3])
    # unitcell[5] = unitcell[1]

    return unitcell
end



# Generates a vector of Suzuki-Trotter gates for a Hubbard Hamiltonian
function gen_gates(sites::Vector{Index{Int64}}, dtau::Float64, operator::String; secondorder = true, U = 0.0, μ = 0.0, finite = false)

    gates = [ITensor() for _ in 1:2] # Define gates as an empty Vector{ITensor}

    if finite
        N += -1
    end
    
    for i in 1:2

        h = get_operator(operator, sites, i; U = U, μ = μ)

        # Second-order ST ordering: dτ is divided by two for odd gates, which will be applied twice
        if secondorder
            if isodd(i)
                gate = exp(-dtau / 2 * h)
            else
                gate = exp(-dtau * h)
            end
        else
            gate = exp(-dtau * h)
        end

        gates[i] = gate 
    end

    return gates
end



# Finds the inverse of a square diagonal two-index tensor
function inv_diagonal(S::ITensor)

    Sinv = copy(S)
    s1, s2 = inds(Sinv)
    values = Int(size(S)[1])
    
    for k in 1:values
        Sinv[s1=>k, s2=>k] = 1 / (S[s1=>k, s2=>k])
    end

    return Sinv
end



# Applies a two-site gate to an iMPS of type Vector{ITensor}, in Vidal form and with a two-site unit cell, and returns the split iMPS sites
function apply_gate(unitcell::Vector{ITensor}, gate::ITensor; odd = true, cutoff = 1e-10, maxdim = 30, plot_coeffs = false, step=1)

    # Pick one of the two sites to time-evolve
    if odd
        site_ind = 2
    else
        site_ind = 4
    end
    
    nextsite_ind = mod1(site_ind + 2, 4)


    # Find physical site indices, for later
    l1, p1, r1 = inds(unitcell[site_ind])
    l2, p2, r2 = inds(unitcell[nextsite_ind])

    # Contract over all iMPS where the central Λ is either Λ_2 (odd gate) or Λ_1 (even gate, considering a Λ_1 - Γ_1 ordering)
    prodsite, lat_clone = contract_full_iMPS(unitcell, central = Int(site_ind/2))

    # Contract the joint site with the two-site operator
    prodsite *= gate


    # Define combiner tensors for left and right sides: this will help switch prodsite_evol
    # from a D,d,d,D indexed tensor to a D*d, D*d tensor, enabling SVD
    i1, i2, i3, i4 = inds(prodsite)

    # ITensor shenanigans
    if !odd
        prodsite = permute(prodsite, i1, i2, i4, i3)
        i1, i2, i3, i4 = inds(prodsite)
    end

    #prodsite = permute(prodsite, i1, i3, i2, i4)
    # println("Inds prodsite: ", inds(prodsite))

    cleft = combiner(i1, i3)
    cright = combiner(i2, i4)
    # println("Inds cleft: ", inds(cleft))
    # println("Inds cright: ", inds(cright))

    # Apply combiners, making prodsite_matrix a D*d, D*d tensor
    prodsite_matrix = prodsite * cleft
    prodsite_matrix *= cright


    # Perform SVD
    U, S, V = svd(prodsite_matrix, inds(prodsite_matrix)[1]; maxdim = maxdim)

    
    # S will be saved as the new central Λ matrix between these two tensors
    # This is where we can normalize the iMPS!
    S /= norm(S)
    S_Lind, S_Rind = inds(S)
    unitcell[mod1(site_ind + 1, 4)] = replaceinds(S, (S_Lind => r1, S_Rind => l2))


    # Find the inverse of the Λ matrix that isn't central in this gate operation to remove it from the sides
    inv_lateral = inv_diagonal(lat_clone)
    println("AQiiiii", norm(inv_lateral - lat_clone))
    latΛ_Lind, latΛ_Rind = inds(inv_lateral)
    inv_lateral_left = replaceind(inv_lateral, latΛ_Lind => l1)
    inv_lateral_right = replaceind(inv_lateral, latΛ_Rind => r2)


    # Recover physical/internal bonds
    U *= dag(cleft) 
    V *= dag(cright)

    U_Lind, U_pind, U_Rind = inds(U)
    V_Lind, V_pind, V_Rind = inds(V)

    # Replace automatic SVD indices by defined internal ones, and 
    replaceinds!(U, (U_Rind => r1, U_Lind => latΛ_Rind)) # Internal bonds have been replaced by SVD bonds, still of χ dimension
    replaceinds!(V, (V_Lind => l2, V_Rind => latΛ_Lind))


    # println("INDS U: ", inds(U))
    # println("INDS V: ", inds(V))

    unitcell[site_ind] = noprime(inv_lateral_left * U)
    unitcell[nextsite_ind] = noprime(V * inv_lateral_right)
    
    # println("INDS LEFT SITE: ", inds(unitcell[site_ind]))
    # println("INDS RIGHT SITE: ", inds(unitcell[nextsite_ind]))

    return unitcell
    
end


# Applies a Suzuki-Trotter step to an iMPS, returns the prodsite matrix to spare the contraction
# for expectation value computation
function ST_step(unitcell::Vector{ITensor}, gates::Vector{ITensor}; secondorder = true, maxdim = 30, plot_coeffs = true)

    ngates = 2

    # Odd sites
    unitcell = apply_gate(unitcell, gates[1], odd = true, maxdim = maxdim)
    
    # Even sites
    unitcell = apply_gate(unitcell, gates[2], odd = false, maxdim = maxdim)
    
    # In second-order ST evolution odd gates are applied again
    if secondorder
        unitcell = apply_gate(unitcell, gates[1], odd = true, maxdim = maxdim)
    end

    return unitcell
end



# Returns the matrix of a given jth-site operator, commonly used in this work
function get_operator(name::String, sites::Vector{Index{Int64}}, site::Int64; t = 1.0, U = 1.0, μ = 0.0)

    N = length(sites)
    nextsite = mod1(site + 1, N) # Index of next site

    if name == "TunnellingUP"
        operator = op("Cdagup", sites[site]) * op("Cup", sites[nextsite])
        operator += op("Cup", sites[site]) * op("Cdagup", sites[nextsite])
        operator *= -t

    elseif name == "TunnellingDOWN"
        operator = op("Cdagdn", sites[site]) * op("Cdn", sites[nextsite])
        operator += op("Cdn", sites[site]) * op("Cdagdn", sites[nextsite])
        operator *= -t

    elseif name == "Tunnelling"
        operator = op("Cdagup", sites[site]) * op("Cup", sites[nextsite])
        operator += op("Cup", sites[site]) * op("Cdagup", sites[nextsite])
        operator += op("Cdagdn", sites[site]) * op("Cdn", sites[nextsite])
        operator += op("Cdn", sites[site]) * op("Cdagdn", sites[nextsite])
        operator *= -t
    
    elseif name == "Onsite"
        operator = U * op("Nupdn", sites[site]) * op("Id", sites[nextsite])

    elseif name == "ChemPot"
        operator = -μ * op("Ntot", sites[site]) * op("Id", sites[nextsite])
    
    elseif name == "Hubbard"
        operator = get_operator("Tunnelling", sites, site)
        #operator2 = U * op("Nupdn", sites[site])
        operator += U * op("Nupdn", sites[site]) * op("Id", sites[nextsite])
    
    elseif name == "GCHubbard"
        operator = get_operator("Tunnelling", sites, site)
        operator += U * op("Nupdn", sites[site]) * op("Id", sites[nextsite])
        operator += -μ * op("Ntot", sites[site]) * op("Id", sites[nextsite])
    
    else
        error("\n\n No operator identified with the following name: ", name)
        #throw(InterruptException())
    end

    return operator
end


# Computes the expectation value of an iMPS for a two-site unit cell in Vidal form
function imps_expect_vidal(unitcell::Vector{ITensor}, operator::ITensor; tolerance = 1e-8)

    # Contract over all tensors in the unit cell: according to this program's construction this accounts for lateral environment effects
    # (Λ_2 is duplicate in the unitcell vector, and it suffices to represent the infinite environment in Vidal form)
    contracted_iMPS, lat_clone = contract_full_iMPS(unitcell, central = 2)
    println("contraaaa", inds(contracted_iMPS))
    leftind, p1, p2, rightind = inds(contracted_iMPS)

    # Absorb the operator into the iMPS contracted cells
    # Fill with identities if necessary to make the operator match all physical indices
    observable = contracted_iMPS * operator

    # Multiply by adjoint mps, ⟨ψ|o|ψ⟩
    observable *= adjoint(contracted_iMPS)

    # Finally contract over lateral indices to obtain scalar
    left_identity = op("Id", leftind)
    right_identity = op("Id", rightind)

    observable *= left_identity
    observable *= right_identity

    # Return observable
    return scalar(observable)
end



# Contracts over all tensors of an iMPS to prepare for multiple observable computation. Allows to choose central Λ in case of iTEBD
function contract_full_iMPS(psi::Vector{ITensor}; central = 1)

    unitcell = copy(psi)
    
    # Find indices of all sites, for comfort
    l1_lamb, r1_lamb = inds(unitcell[1])
    l1, p1, r1 = inds(unitcell[2])
    l2_lamb, r2_lamb = inds(unitcell[3])
    l2, p2, r2 = inds(unitcell[4])

    if central == 1
        lateral_clone = unitcell[3]

        # Replace indices to allow for different iMPS contraction
        unitcell[3] = prime(unitcell[3], l2_lamb) # Avoid looping, now \Lambda_2 is the tensor on both sides
        unitcell[4] = replaceinds(unitcell[4], (r2 => l1_lamb))

        # Output tensor will have internal indices of l2_lamb', r2_lamb
        # Output lateral_clone will have indices r1, l2_lamb
    else
        lateral_clone = unitcell[1]
        
        # Prepare the lateral clone to be absorbed through the left
        lateral_clone = replaceinds(lateral_clone, (l1_lamb => r2))
        unitcell[1] = prime(unitcell[1], l1_lamb) # For consistency with previous result

        # Output tensor will have internal indices of l1_lamb', r1_lamb
        # Output lateral_clone will have indices r2, l1_lamb
    end
 

    # Contract over all tensors in the unit cell: according to this program's construction this accounts for lateral environment effects
    # (Λ_2 is duplicate in the unitcell vector, and it suffices to represent the infinite environment in Vidal form)

    # Carefully picking the order of contraction will output a tensor with the most comfortable index ordering
    # (virtual-physical-physical-virtual)
    init_contract = mod(central, 2)*2 + 1
    contracted_iMPS = unitcell[init_contract]

    for i in 1:3
        nextsite = mod1(init_contract + i, 4)
        contracted_iMPS *= unitcell[nextsite]
    end

    # Finally contract over the lateral clone at the rightmost position
    contracted_iMPS *= lateral_clone

    println("INDS CONTRACTED HERE !!! ", inds(contracted_iMPS))
    

    # Returns a fully internally-contracted iMPS unit cell
    return contracted_iMPS, lateral_clone
end





# Computes the expectation value of an iMPS for a two-site unit cell in Vidal form. Spares tensor contraction in case of more than one observable
function imps_expect_vidal_optimal(contracted_iMPS::ITensor, operator::ITensor)

    leftind, p1, p2, rightind = inds(contracted_iMPS)

    # Absorb the operator into the iMPS contracted cells
    # Fill with identities if necessary to make the operator match all physical indices
    observable = contracted_iMPS * operator

    # Multiply by adjoint mps, ⟨ψ|o|ψ⟩
    observable *= adj(contracted_iMPS)

    # Finally contract over lateral indices to obtain scalar
    left_identity = op("Id", leftind)
    right_identity = op("Id", rightind)

    observable *= left_identity
    observable *= right_identity

    # Return observable
    return scalar(observable)
end




function loop_iMPS(unitcell::Vector{ITensor})

    N = length(unitcell)

    # Get lateral indices
    l1, p1, r1 = inds(unitcell[1])
    ln, pn, rn = inds(unitcell[N])

    replaceinds!(unitcell[1], (l1 => rn))

    return unitcell
end


function unloop_iMPS(unitcell::Vector{ITensor}, leftindex::Index{Int64})

    # Get first cell indices
    rn, p1, r1 = inds(unitcell[1])

    replaceinds!(unitcell[1], (rn => leftindex))

    return unitcell
end



# Computes the Von Neumann entropy of an iMPS at a given site
function site_entropy(iMPS::Vector{ITensor}, site::Int64)

    

    entropy = 0


    return entropy
end





let 
    # 
    bdim = 16

    # Create an iMPS with a known initial state
    psi, sites, links = init_iMPS([1,1], bdim = bdim)

    # Define hamiltonian parameters
    U = 1.0
    μ = 0.0
    operator = "Hubbard"

    # Simulation parameters
    #   iTEBD
    cutoff = 1e-5
    dtau = 0.01
    steps = 2000
    #   Result analysis
    checkevery = 20
    framespersecond = 10
    tolerance = 1e-7
    mirar_correlacions = false
    finite = false
    secondorder = true

    mesures = []
    lind = inds(psi[1])[1]

    # Generate iTEBD gates
    gates = gen_gates(sites, dtau, operator, U = U, μ = μ, secondorder = secondorder)

    # Prepare iMPS to apply iTEBD
    psi = normalise_unitcell(psi)
    #psi = loop_iMPS(psi)

    # probability vector, should be 1 through all iterations for particle number conservation
    totalupprob = []
    totaldnprob = []

    #psi = ST_step(psi, gates, maxdim = bdim, secondorder = secondorder)
    #psi = apply_gate(psi, gates[1], odd = false, maxdim = bdim)

    println("INDEXOS DE PSI")
    println("Tensor 1: ", inds(psi[1]))
    println("Tensor 2: ", inds(psi[2]))
    println("Tensor 3: ", inds(psi[3]))
    println("Tensor 4: ", inds(psi[4]))


    op1 = op("Id", sites[1]) * op("Id", sites[2])
    prob1 = imps_expect_vidal(psi, op1)

    println("Ara mirem escalar: ", inds(prob1))

    # Begin iTEBD loop
    # for step in 1:steps

        # Evolve the system using iTEBD, second-order Suzuki-Trotter gate ordering:
        # First (square-rooted) odd gates are applied, then even, then odd gates again
        # The functions gen_gates and secondorder_STstep defined in the file "functions.jl"
        # implement this automatically, generating an array in which the gates are well-ordered
        #psi = ST_step(psi, gatetest, maxdim = bdim, cutoff=cutoff, step = step, secondorder = secondorder, plot_coeffs = mirar_correlacions, finite = finite)


        #if mod(step, checkevery) == 0

        #    println("Step ", step, " of ", steps)
            
            # Measure density profile
        #    dens_up = zeros(N)
        #    dens_dn = zeros(N)

            # Contract over all sites to compute observable
        #    contracted_iMPS = contract_full_iMPS(psi)


        #    for pos in 1:2
        #        op_dens_up = op("Nup", sites[pos])
        #        op_dens_down = op("Ndn", sites[pos])

#                othersite = mod1(pos + 1, 2)
#                op_dens_up *= op("Id", sites[othersite])
#                op_dens_down *= op("Id", sites[othersite])

                # Find expected values
#                dens_up[pos] = abs(imps_expect_vidal_optimal(contracted_iMPS, op_dens_up))
#                dens_dn[pos] = abs(imps_expect_vidal_optimal(contracted_iMPS, op_dens_dn))
#            end

#            push!(totalupprob, sum(dens_up))
#            push!(totaldnprob, sum(dens_dn))

#            plot(dens_up, ylims=(0,1), legend=false, title="* Iteració = $step", lc=:red, linewidth=3) # Nup
#            plot!(dens_dn, ylims=(0,1), legend=false, lc=:blue, linewidth=3) # Ndn
#            frame(anim)

#        end


    #end




    #j = 1
    #gate1 = op("Ntot", sites[j])
    #expect2 = scalar(imps_expect(psi, gate1, j))
    #expect2 = imps_expect(psi, gate1, 2)
    #print("Escalar final = ", expect2, "\n")
    #print("Escalar final = ", inds(expect2), "\n")
    #print()

    # Observable measurement
    #n = expect(psi, "Nup")

    # Generate density plot
    #plot(n, ylimits = (0, 1))
    #savefig("density_1part.png")

    # Generate iTEBD animation
#    prelude = "figures_vidal/finite_"*string(finite)*"/N_"*string(N)*"/op_"*operator*"/"
#    fig_id = "U_"*string(U)*"_steps"*string(steps)
#    gif(anim, prelude*"itebdanim"*fig_id*".gif", fps=framespersecond)


#    plot(totalupprob, ylim=(0,3), lc=:red, title="# de partícules")
#    plot!(totaldnprob, ylim=(0,3), lc=:blue, title="# de partícules")
#    savefig(prelude*"Npart"*fig_id*".png")

#    if mirar_correlacions
#        gif(anim2, prelude*"Correlacions"*fig_id*".gif", fps=framespersecond*5)
#    end

#    println("Hamiltonian: ", operator, " | N = ", N, " | finite = ", finite)

    # Generate energy evolution plot
    #plot(energies, yformatter = :scientific, ylimits=(-1e-15, 1e-15))
    #savefig("energy_evol_itebd.png")
 
end