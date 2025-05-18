# iTEBD evolution using vidal notation for the unit cell
using ITensors, ITensorMPS
using Plots
gr()

#include("functions.jl")
anim = Animation()
anim2 = Animation()

# Generates an iMPS ITensor in Vidal form with a given initial state. 1:empty, 2:up, 3:dn, 4:updn
function init_iMPS(state::Vector; bdim=16)
    
    references = Dict(0=>"Emp", 1=>"Up", 2=>"Dn", 3=>"UpDn")
    
    # Initializes iMPS with two-site unit cell 
    sites = siteinds("Electron", 2)
    links = [Index(bdim, "link-$i") for i in 1:4]
    psi = [ITensor() for _ in 1:4]
    println(sites[1])

    for i in 1:4
        if iseven(i) # Γ
            T = ITensor(links[i], sites[Int(i/2)], links[mod1(i+1, 4)])
            T[links[i]=>1, sites[Int(i/2)]=>references[state[Int(i/2)]], links[mod1(i+1, 4)]=>1] = 1.0
            psi[i] = T
            #psi[i] = randomITensor(links[i], sites[Int(i/2)], links[i+1])

        else # Λ
            psi[i] = ITensor(zeros(bdim, bdim), links[i], links[mod1(i+1, 4)])
            for k in 1:bdim
                psi[i][links[i]=>k, links[mod1(i+1, 4)]=>k] = 1.0
            end
        end
    end


    k = 3
    #psi[1][links[1]=>k, links[2]=>k]=1.0
    psi = normalise_unitcell(psi, bdim = bdim)
    println("Norm of diag tensor: ", sumdiag(psi[1]))

    return psi, sites, links
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


# Contracts over all tensors of an iMPS to prepare for multiple observable computation. Allows to choose central Λ in case of iTEBD
function contract_full_iMPS(psi::Vector{ITensor}; central = 1)

    unitcell = copy(psi)
    # New attempt: the link indices are defined in a cyclic way, spares contraction conditions

    # Central == 1 means that the matrix Λ_1 is taken to be the centre of orthogonality of this contraction

    if central == 1     # Λ_2 - Γ_2 - Λ_1 - Γ_1 - Λ_2
        rightmost_clone = copy(unitcell[3])

        unitcell[3] = prime(unitcell[3], commonind(unitcell[2], unitcell[3])) # Leftmost contraction
        rightmost_clone = prime(rightmost_clone, commonind(rightmost_clone, unitcell[4]))

        # Output tensor will have internal indices of l2_lamb', r2_lamb'

    else                # Λ_1 - Γ_1 - Λ_2 - Γ_2 - Λ_1
        rightmost_clone = copy(unitcell[1])
        
        unitcell[1] = prime(unitcell[1], commonind(unitcell[1], unitcell[4]))
        rightmost_clone = prime(rightmost_clone, commonind(rightmost_clone, unitcell[2]))

        # Output tensor will have internal indices of l1_lamb', r1_lamb'
    end
 

    # Contract over all tensors in the unit cell: according to this program's construction this accounts for lateral environment effects
    # (Λ_2 is duplicate in the unitcell vector, and it suffices to represent the infinite environment in Vidal form)
    contracted_iMPS = unitcell[1]

    for i in 2:4
        contracted_iMPS *= unitcell[i]
    end

    # Finally contract over the lateral clone at the rightmost position
    contracted_iMPS *= rightmost_clone

    # Returns a fully internally-contracted iMPS unit cell
    return contracted_iMPS, rightmost_clone
end

function normalise_unitcell(unitcell::Vector{ITensor}; bdim = 16)

    #unitcell[1] = normalise_diagonal(unitcell[1])
    #unitcell[3] = normalise_diagonal(unitcell[3])
    unitcell[1] /= norm(unitcell[1])
    unitcell[3] /= norm(unitcell[3])
    #unitcell[2] /= norm(unitcell[2])
    #unitcell[4] /= norm(unitcell[4])
    # unitcell[5] = unitcell[1]
    factor = 4*sqrt(2)*2^(1/6)#2#^(5/6)*2^(1/6)
    unitcell[1] *= factor
    unitcell[3] *= factor

    return unitcell
end


function sumdiag(S::ITensor)

    mida = Int(size(S)[1])
    S_Lind, S_Rind = inds(S)

    sumatot = 0.0
    for k in 1:mida
        sumatot += S[S_Lind=>k, S_Rind=>k]^2
    end
    return sumatot
end


# Normalises a square diagonal two-index tensor
function normalise_diagonal(S::ITensor)

    mida = Int(size(S)[1])
    S_Lind, S_Rind = inds(S)

    sumatot = 0.0
    for k in 1:mida
        sumatot += S[S_Lind=>k, S_Rind=>k]^2
    end

    return S / sqrt(sumatot)
end

# Finds the site index of a Γ tensor
function find_siteindex(sitetensor::ITensor)

    for index in inds(sitetensor)
        if hastags(index, "Site")
            return index
        end
    end
    return
end

# Converts a prodsite tensor (D,d,d,D) to matrix for SVD
function tensor_to_matrix(unitcell::Vector{ITensor}, rightmost_clone::ITensor, central::Int64, prodsite::ITensor)

    if central == 1
        #left_internal = prime(commonind(unitcell[2], rightmost_clone)) # Primed, as we have just done
        left_internal = commonind(unitcell[2], rightmost_clone)
        left_physical = find_siteindex(unitcell[4])
        #right_internal = prime(commonind(unitcell[3], unitcell[4]))
        right_internal = commonind(unitcell[3], unitcell[2])
        right_physical = find_siteindex(unitcell[2])

        cleft = combiner(left_internal, left_physical)
        cright = combiner(right_internal, right_physical)

    else
        #left_internal = prime(commonind(unitcell[4], rightmost_clone)) # Primed, as we have just done
        left_internal = commonind(unitcell[4], rightmost_clone)
        left_physical = find_siteindex(unitcell[2])
        #right_internal = prime(commonind(unitcell[1], unitcell[2]))
        right_internal = commonind(unitcell[1], unitcell[2])
        right_physical = find_siteindex(unitcell[4])

        cleft = combiner(left_internal, left_physical)
        cright = combiner(right_internal, right_physical)
    end

    return prodsite * cleft * cright, cleft, cright
end


# Generates a vector of Suzuki-Trotter gates for a Hubbard Hamiltonian
function gen_gates(sites::Vector{Index{Int64}}, dtau::Float64, operator::String; secondorder = true, U = 0.0, μ = 0.0, V = 0.0, finite = false)

    gates = [ITensor() for _ in 1:2] # Define gates as an empty Vector{ITensor}
    
    for i in 1:2

        h = get_operator(operator, sites, i; U = U, μ = μ, V = V)

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


# Applies a two-site gate to an iMPS of type Vector{ITensor}, in Vidal form and with a two-site unit cell, and returns the split iMPS sites
function apply_gate(unitcell::Vector{ITensor}, gate::ITensor; odd = true, cutoff = 1e-10, maxdim = 30, plot_coeffs = false, step=1)

    # Pick one of the two sites to time-evolve
    odd = true
    if odd
        site_ind = 2 # In the odd gate application the sites are contracted together over Λ_2
    else
        site_ind = 4
    end
    nextsite_ind = mod1(site_ind + 2, 4)

    # Contract unit cell onto a single tensor of dimension D,d,d,D
    prodsite, rightmost_clone = contract_full_iMPS(unitcell, central = mod1(site_ind+1,4))

    # Apply Trotter gate to contracted unit cell
    prodsite *= gate
    noprime!(prodsite)

    # Convert product tensor to matrix of dimension D*d, D*d
    prodsite_matrix, cleft, cright = tensor_to_matrix(unitcell, rightmost_clone, mod1(site_ind+1,4), prodsite)

    # Perform SVD
    U, S, V = svd(prodsite_matrix, commonind(prodsite_matrix, cleft), maxdim = maxdim, cutoff = cutoff)

    # S will be saved as the new central Λ matrix between these two tensors
    # This is where we can normalize the iMPS!
    S = normalise_diagonal(S)

    # SVD indices are conserved to allow for variable internal bond dimension
    unitcell[mod1(site_ind + 1, 4)] = S

    # Recover physical/internal bonds
    U *= dag(cleft)
    noprime!(U)
    V *= dag(cright)
    noprime!(V)

    # Find the inverse of the Λ matrix that isn't central in this gate operation to remove it from the sides
    # Remove primed indices to enable contraction
    noprime!(rightmost_clone)
    inv_lateral = inv_diagonal(rightmost_clone)
    inv_lateral = normalise_diagonal(inv_lateral)
    
    # Update Γ tensors
    unitcell[site_ind] = U * inv_lateral
    unitcell[nextsite_ind] = V * inv_lateral

    return unitcell
end


# Applies a Suzuki-Trotter step to an iMPS, returns the prodsite matrix to spare the contraction
# for expectation value computation
function ST_step(unitcell::Vector{ITensor}, gates::Vector{ITensor}; secondorder = true, maxdim = 30, cutoff = 1e-10, plot_coeffs = true, bdim = 16)

    ngates = 2

    # Odd sites
    unitcell = apply_gate(unitcell, gates[1], odd = true, maxdim = maxdim, cutoff = cutoff)
    unitcell = normalise_unitcell(unitcell, bdim = bdim)
    

    #println("Inds 1: ", inds(unitcell[1]))
    #println("Inds 2: ", inds(unitcell[2]))
    #println("Inds 3: ", inds(unitcell[3]))
    #println("Inds 4: ", inds(unitcell[4]))
    
    # Even sites
    unitcell = apply_gate(unitcell, gates[2], odd = false, maxdim = maxdim, cutoff = cutoff)
    unitcell = normalise_unitcell(unitcell, bdim = bdim)
    
    # In second-order ST evolution odd gates are applied again
    if secondorder
        unitcell = apply_gate(unitcell, gates[1], odd = true, maxdim = maxdim, cutoff = cutoff)
        unitcell = normalise_unitcell(unitcell, bdim = bdim)
    end

    return unitcell
end



# Returns the matrix of a given jth-site operator, commonly used in this work
function get_operator(name::String, sites::Vector{Index{Int64}}, site::Int64; t = 1.0, U = 1.0, μ = 0.0, V = 0.0)

    N = 2
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
        operator = -t * op("Cdagup", sites[site]) * op("Cup", sites[nextsite])
        operator += -t * op("Cup", sites[site]) * op("Cdagup", sites[nextsite])
        operator += -t * op("Cdagdn", sites[site]) * op("Cdn", sites[nextsite])
        operator += -t * op("Cdn", sites[site]) * op("Cdagdn", sites[nextsite])
        #operator *= -t
    
    elseif name == "Onsite"
        operator = U * op("Nupdn", sites[site]) * op("Id", sites[nextsite])

    elseif name == "ChemPot"
        operator = -μ * op("Ntot", sites[site]) * op("Id", sites[nextsite])
    
    elseif name == "Hubbard"
        operator = get_operator("Tunnelling", sites, site)
        #operator2 = U * op("Nupdn", sites[site])
        operator += U * op("Nupdn", sites[site]) * op("Id", sites[nextsite])
        operator += U * op("Id", sites[site]) * op("Nupdn", sites[nextsite])
    
    elseif name == "GCHubbard"
        operator = get_operator("Tunnelling", sites, site)
        operator += U * op("Nupdn", sites[site]) * op("Id", sites[nextsite])
        operator += U * op("Id", sites[site]) * op("Nupdn", sites[nextsite])
        operator += -μ * op("Ntot", sites[site]) * op("Id", sites[nextsite])
        operator += -μ * op("Id", sites[site]) * op("Ntot", sites[nextsite])
    
    elseif name == "ExtHubbard"
        operator = get_operator("Tunnelling", sites, site)
        operator += U * op("Nupdn", sites[site]) * op("Id", sites[nextsite])
        operator += U * op("Id", sites[site]) * op("Nupdn", sites[nextsite])
        operator += V * op("Ntot", sites[site]) * op("Ntot", sites[nextsite])
    
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
    contracted_iMPS, lat_clone = contract_full_iMPS(unitcell, central = 1)

    # Remove primed indices
    noprime!(contracted_iMPS)

    # Absorb the operator into the iMPS contracted cells
    # Fill with identities if necessary to make the operator match all physical indices
    observable = contracted_iMPS * operator

    # Multiply by adjoint mps, ⟨ψ|o|ψ⟩
    observable *= adjoint(contracted_iMPS)

    # Finally contract over lateral indices to obtain scalar
    leftind, rightind = commoninds(observable, contracted_iMPS)
    left_identity = op("Id", leftind)
    right_identity = op("Id", rightind)

    observable *= left_identity
    observable *= right_identity

    # Return observable
    return scalar(observable)
end



# Computes the expectation value of an iMPS for a two-site unit cell in Vidal form. Spares tensor contraction in case of more than one observable
function imps_expect_vidal_optimal(contracted_iMPS::ITensor, operator::ITensor)

    # Remove primed indices
    noprime!(contracted_iMPS)

    # Absorb the operator into the iMPS contracted cells
    # Fill with identities if necessary to make the operator match all physical indices
    observable = contracted_iMPS * operator

    # Multiply by adjoint mps, ⟨ψ|o|ψ⟩
    observable *= adjoint(contracted_iMPS)

    # Finally contract over lateral indices to obtain scalar
    leftind, rightind = commoninds(observable, contracted_iMPS)
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
    # N=2
    bdim = 16
    println("Bond dimension: ", bdim)

    # Create an iMPS with a known initial state
    psi, sites, links = init_iMPS([1,0], bdim = bdim)

    # Define hamiltonian parameters
    U = 10.0
    μ = 0.0
    V = 5.0
    operator = "Hubbard"

    # Simulation parameters
    #   iTEBD
    cutoff = 1e-10
    dtau = 0.01
    steps = 1500
    #   Result analysis
    checkevery = 50
    framespersecond = 10
    tolerance = 1e-7
    mirar_correlacions = false
    finite = false
    secondorder = true
    savedata = true

    mesures = []
    lind = inds(psi[1])[1]
    #println(psi[1])

    # Generate iTEBD gates
    gates = gen_gates(sites, dtau, operator, U = U, μ = μ, V = V, secondorder = secondorder)
    

    # Prepare iMPS to apply iTEBD
    #psi = normalise_unitcell(psi)

    # probability vector, should be 1 through all iterations for particle number conservation
    totalupprob = []
    totaldnprob = []

    # psi = apply_gate(psi, gates[1], odd = true, maxdim = bdim, cutoff = cutoff)
    #  println(psi[2])
    psi = normalise_unitcell(psi, bdim = bdim)

    op_meas_energy = get_operator(operator, sites, 1, U=U, μ = μ, V = V)
    upsite = op("Nup", sites[1]) * op("Id", sites[2])
    dnsite = op("Ndn", sites[2]) * op("Id", sites[1])

    prob1 = imps_expect_vidal(psi, upsite)
    prob2 = imps_expect_vidal(psi, dnsite)
    ene1 = imps_expect_vidal(psi, op_meas_energy)

    println("Ara mirem escalar: ", prob1)
    println("Ara mirem escalar: ", prob2)
    println("Ara mirem energia: ", ene1)
    #println("Ara mirem escalar: ", prob2)


    j = 1
    #op_meas_energy = get_operator("GCHubbard", sites, 1, U=U, μ = μ, V = V)
    energies = []


    # Begin iTEBD loop
    for step in 1:steps

        # Evolve the system using iTEBD, second-order Suzuki-Trotter gate ordering:
        # First (square-rooted) odd gates are applied, then even, then odd gates again
        # The functions gen_gates and secondorder_STstep defined in the file "functions.jl"
        # implement this automatically, generating an array in which the gates are well-ordered
        psi = ST_step(psi, gates, maxdim = bdim, cutoff = cutoff, secondorder = secondorder, bdim = bdim)
        psi = normalise_unitcell(psi, bdim = bdim)

        if mod(step, checkevery) == 0

            println("Step ", step, " of ", steps)
            
            # Measure density profile
            dens_up = zeros(2)
            dens_dn = zeros(2)
            
            contracted_iMPS, lateral_clone = contract_full_iMPS(psi, central = 2)

            energy = imps_expect_vidal_optimal(contracted_iMPS, op_meas_energy)


            for pos in 1:2
                # Contract over all sites to compute observable
                op_dens_up = op("Nup", sites[pos])
                op_dens_down = op("Ndn", sites[pos])

                othersite = mod1(pos + 1, 2)
                op_dens_up *= op("Id", sites[othersite])
                op_dens_down *= op("Id", sites[othersite])

                # Find expected values
                dens_up[pos] = imps_expect_vidal_optimal(contracted_iMPS, op_dens_up)
                dens_dn[pos] = imps_expect_vidal_optimal(contracted_iMPS, op_dens_down)
            end

            push!(energies, energy)
            push!(totalupprob, sum(dens_up))
            push!(totaldnprob, sum(dens_dn))

            plot(dens_up, ylims=(0,1), legend=false, title="* Iteració = $step", lc=:red, linewidth=3) # Nup
            plot!(dens_dn, ylims=(0,1), legend=false, lc=:blue, linewidth=3) # Ndn
            frame(anim)

        end


    end




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
    if savedata
        prelude = "figures_vidal/op_"*operator*"/"
        fig_id = "U_"*string(U)*"_V_"*string(V)*"_dt_"*string(dtau)*"_bdim_"*string(bdim)*"_steps"*string(steps)
        gif(anim, prelude*"itebdanim"*fig_id*".gif", fps=framespersecond)


        plot(totalupprob, ylim=(0,2), lc=:red, title="# de partícules")
        plot!(totaldnprob, ylim=(0,2), lc=:blue, title="# de partícules")
        savefig(prelude*"Npart"*fig_id*".png")


        plot(energies, lc=:blue, title="Energia")
        savefig(prelude*"Energy"*fig_id*".png")

        #    if mirar_correlacions
        #        gif(anim2, prelude*"Correlacions"*fig_id*".gif", fps=framespersecond*5)
        #    end
    end

    println("Hamiltonian: ", operator, " | finite = ", finite)

    # Generate energy evolution plot
    #plot(energies, yformatter = :scientific, ylimits=(-1e-15, 1e-15))
    #savefig("energy_evol_itebd.png")
 
end