using ITensors, ITensorMPS
using Plots
gr()

#include("functions.jl")
anim = Animation()



# Generates a FH Hamiltonian MPO for U, μ values
function FH_Hamiltonian(N, sites, U, μ)
    t = 1.0 # Hopping amplitude
    os = OpSum()
    for j in 1:(N - 1)

        # Hopping terms
        os += -t, "Cdagup", j, "Cup", mod1(j + 1, N)
        os += -t, "Cup", j, "Cdagup", mod1(j + 1, N)
        os += -t, "Cdagdn", j, "Cdn", mod1(j + 1, N)
        os += -t, "Cdn", j, "Cdagdn", mod1(j + 1, N)

        # On-site interaction
        os += U, "Nup", j, "Ndn", j

        # Chemical potential
        os += μ, "Ntot", j
    end

    os += U, "Nup", N, "Ndn", N
    os += μ, "Ntot", N

    H = MPO(os, sites)

    return H
end


# Initializes an MPS as a Vector{ITensor}. Allows for iMPS selection
function init_iMPS(N::Int64, sites::Vector, links::Vector; infinite = true)

    psi = [randomITensor(links[mod1(i-1, N)], sites[i], links[i]) for i in 1:N]

    if !infinite
        psi[1] = randomITensor(sites[1], links[1])
        psi[N] = randomITensor(links[N-1], sites[N])
    end

    return psi
end


# Sets a given FH state for an MPS. Input is the MPS and a vector where 1:empty, 2:up, 3:dn, 4:updn
function set_FHstate(iMPS::Vector{ITensor}, sites::Vector, links::Vector, state::Vector; bdim=10)

    N = length(iMPS)

    references = Dict(0=>"Emp", 1=>"Up", 2=>"Dn", 3=>"UpDn")

    for i in 1:N
        #print("AAA", typeof(state[i]), "\n")
        T = ITensor(links[i], sites[i], links[mod1(i+1,N)])
        T[links[i]=>bdim, sites[i]=>references[state[i]], links[mod1(i+1,N)]=>bdim] = 1.0
        iMPS[i] = T
        #MPS[i][links(i), sites[i]=>Dict(state[i]), links(mod1(i+1,N))] = 1.0
    end

    return iMPS
end


# ST ORDERING, EXPLAIN
function ST_indexing(pos::Int64, N::Int64)

    if mod(N,2) == 1
        ST_index = mod1(2*pos - 1, N)

    # If pos is in an even site, considering ST starts in an odd site:
    elseif pos > Int64(N/2)
        ST_index = mod1(2*pos, N)
    else
        ST_index = mod1(1 + 2*(pos - 1), N)
    end

    return ST_index
end


# Generates a vector of Suzuki-Trotter gates for a Hubbard Hamiltonian
function gen_gates(N::Int64, sites::Vector{Index{Int64}}, dtau::Float64; secondorder = True)

    #gates = [ITensor() for _ in 1:N] # Define gates as an empty Vector{ITensor}
    gates = []

    for i in 1:N

        # Periodic boundary conditions
        iplusone = mod1(i+1, N)

        h = get_operator("Hubbard", sites, i)

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

        #if isodd(i)
        #    push!(odd_gates, gate)
        #else
        #    push!(even_gates, gate)
        #end

        push!(gates, gate)

        # gates are saved in the gates vector in the S-T order, so that a simple iterated application
        # of gates over the iMPS gives already Suzuki-Trotter ordering
        #gatepos = ST_index(N, i)
        #gates[gatepos] = gate
        
    end

    #return odd_gates, even_gates
    return gates
end


# Applies a two-site gate to an iMPS of type Vector{ITensor}, and returns the split iMPS sites
function apply_gate(site1::ITensor, site2::ITensor, gate::ITensor; cutoff = 1e-20, maxdim = 100)

    N = length(iMPS)

    # Contract the site at pos with the next one, two-site gate
    prodsite = site1 * site2

    # Contract the joint site with the two-site operator
    prodsite_evol = gate * prodsite

    # Define combiner tensors for left and right sides: this will help switch prodsite_evol
    # from a D,d,d,D indexed tensor to a D*d, D*d tensor, enabling SVD
    cleft = combiner(inds(prodsite_evol)[1], inds(prodsite_evol)[2])
    cright = combiner(inds(prodsite_evol)[3], inds(prodsite_evol)[4])

    # Apply combiners, making prodsite_matrix a D*d, D*d tensor
    prodsite_matrix = prodsite * cleft
    prodsite_matrix *= cright

    # Perform SVD
    U, S, V = svd(prodsite_matrix, inds(prodsite_matrix)[1])#; maxdim = maxdim)

    # diagonal_array = [S[inds(S)[1]=>i, inds(S)[2]=>i] for i in 1:size(S)[1]]
    
    # Compute square root matrix
    sqS = copy(S)
    for k in 1:size(S)[1]
        sqS[inds(sqS)[1]=>k, inds(sqS)[2]=>k] = sqrt( S[inds(S)[1]=>k, inds(S)[2]=>k] )
    end

    # We have separated S into its two factors, (√S)^2=S: it's a simple operation that can be
    # performed element-wise since we know that S is diagonal and its elements are real.
    # However the two √S matrices that will be absorbed into U and S respectively need to be
    # explicitly connected by a shared index, which will be the link between site1 and site2
    sql = copy(sqS)
    sqr = copy(sqS)
    # Replace the sql right index with the site1-site2 internal index
    replaceindex!(sql, inds(sql)[2], inds(site2)[1])
    replaceindex!(sqr, inds(sqr)[1], inds(site2)[1])
    
    # Finally form and reshape the new site tensors, of indices D,d,D
    L = U * sqSl * dag(cleft)
    R = sqSr * V * dag(cright)

    # Return updated sites
    return L, R
end


# Applies a Suzuki-Trotter step to an MPS
function ST_step(iMPS::Vector{ITensor}, gates::Vector{ITensor}; secondorder = True)

    N = length(iMPS)

    # Odd gates
    for i in 1:N
        if isodd(i)
            iMPS[i], iMPS[mod1(i+1, N)] = apply_gate(iMPS[i], iMPS[mod1(i+1, N)], gates[i])
        end
    end
    
    # Even gates
    for i in 1:N
        if iseven(i)
            iMPS[i], iMPS[mod1(i+1, N)] = apply_gate(iMPS[i], iMPS[mod1(i+1, N)], gates[i])
        end
    end
    
    # In second-order ST evolution odd gates are applied again
    if secondorder
        for i in 1:N
            if isodd(i)
                iMPS[i], iMPS[mod1(i+1, N)] = apply_gate(iMPS[i], iMPS[mod1(i+1, N)], gates[i])
            end
        end
    end

    return iMPS
end


# Returns the matrix of a given jth-site operator, commonly used in this work
function get_operator(name::String, sites::Vector{Index{Int64}}, site::Int64; t = 1.0, U = 0.0, μ = 0.0)

    N = length(sites)
    nsite = mod1(site + 1, N) # Index of next site

    if name == "TunnellingUP"
        operator = op("Cdagup", sites[site]) * op("Cup", sites[mod1(site + 1, N)])
        operator += op("Cup", sites[site]) * op("Cdagup", sites[mod1(site + 1, N)])
        operator *= -t

    elseif name == "TunnellingDOWN"
        operator = op("Cdagdn", sites[site]) * op("Cdn", sites[mod1(site + 1, N)])
        operator += op("Cdn", sites[site]) * op("Cdagdn", sites[mod1(site + 1, N)])
        operator *= -t

    elseif name == "Tunnelling"
        operator = op("Cdagup", sites[site]) * op("Cup", sites[mod1(site + 1, N)])
        operator += op("Cup", sites[site]) * op("Cdagup", sites[mod1(site + 1, N)])
        operator += op("Cdagdn", sites[site]) * op("Cdn", sites[mod1(site + 1, N)])
        operator += op("Cdn", sites[site]) * op("Cdagdn", sites[mod1(site + 1, N)])
        operator *= -t
    
    elseif name == "Onsite"
        operator = -U * op("Nupdn", sites[site])

    elseif name == "ChemPot"
        operator = -μ * op("Ntot", sites[site])
    
    elseif name == "Hubbard"
        operator = op("Cdagup", sites[site])*op("Cup", sites[mod1(site + 1, N)])
        operator += op("Cup", sites[site])*op("Cdagup", sites[mod1(site + 1, N)])
        operator += op("Cdagdn", sites[site])*op("Cdn", sites[mod1(site + 1, N)])
        operator += op("Cdn", sites[site])*op("Cdagdn", sites[mod1(site + 1, N)])
        operator *= -t
        operator += -U * op("Nupdn", sites[site])
    
    else
        print("\n\n No operator identified with the following name: ", name, ". Aborting program \n")
        throw(InterruptException())
    end

    return operator
end


# Returns the expected value of a given operator for an iMPS made of Vector{ITensor}
# Works for operators of any size
# Computes the expected value of an operator in an iMPS Vector{ITensor} object
function imps_expect(iMPS::Vector{ITensor}, operator::ITensor, site::Int64)

    # Find number of sites the operator acts on
    s = Int64(length(size(operator))/2)
    print("\n\t> S = ", s, "\n")

    N = length(iMPS)

    site_tensor = iMPS[site]

    # If the operator is of more than one site, absorb the (s-1) following sites into 
    # a single tensor that can be multiplied by the operator
    if s > 1
        for i in 1:(s-1)
            site_tensor *= iMPS[mod1(site + i, N)]
        end
    end


    # Apply operator gate to sites of interest to create a tensor of dimension D,D,D,D
    applied_tosite = site_tensor' * operator * site_tensor

    # Now a new Vector{ITensor} object is created of length N-(s-1), where s
    # is the number of sites the operator acts upon. In one of its positions
    # the object calculated above will be saved, and in the others the contractions
    # of the sites at both MPS sites (site' * site) will be stored. Finally the
    # whole iMPS will be contracted vertically
    cont_iMPS = [ITensor() for _ in 1:(N - (s - 1))]

    # We save the applied_tosite tensor to the first position of cont_iMPS
    # We don't need to preserve the position of the original iMPS, just the order.
    # cont_iMPS will store first the applied_tosite tensor and then the following
    # iMPS symmetrical sites, contracted, ordered by looping over all the iMPS from the
    # site immediately below the contracted operator up to the one immediately above the
    # site entered as an argument in this function
    cont_iMPS[1] = applied_tosite

    # For each position of this new vector that is not storing the tensor <ϕi|O|ϕi>:
    for item in 1:N-s
        # Set the cursor at position pos: the site acted upon + 1 + (s-1) (sites absorbed by multisite operator)
        pos = mod1(site + item + (s-1), N)
        pos_cont = mod1(item + 1, N - (s-1))
        indsite = "Electron,Site,n="*string(pos) # Physical index tag at given site to unprime

        # Contract symmetric iMPS sites upon which no operator acted and store at pos
        cont_iMPS[pos_cont] = noprime(iMPS[pos]', indsite) * iMPS[pos]
    end

    # Contract rest of iMPS
    expect = cont_iMPS[1] # Any position is valid: all indices are connected and will be contracted
    for i in 2:(N - (s-1))
        expect *= cont_iMPS[i]
    end

    return expect
end


# Computes the Von Neumann entropy of an iMPS at a given site
function site_entropy(iMPS::Vector{ITensor}, site::Int64)

    

    entropy = 0


    return entropy
end




let 
    # 
    N = 5
    bdim = 20
    sites = siteinds("Electron", N)
    links = [Index(bdim, "link-$i") for i in 1:N]

    # Create an iMPS with a known initial state
    psi = init_iMPS(N, sites, links)
    print(typeof(sites), "\n")
    psi = set_FHstate(psi, sites, links, [1, 1, 0, 0, 0])

    # Define simulation parameters
    U = 0.0
    μ = 0.0
    
    # iTEBD parameters, gates and operator
    cutoff = 1e-5
    dtau = 0.01
    steps = 200

    mesures = []

    j = 5

    for step in 1:steps
        
        # Evolve the system using iTEBD, second-order Suzuki-Trotter gate ordering:
        # First (square-rooted) odd gates are applied, then even, then odd gates again
        # The functions gen_gates and secondorder_STstep defined in the file "functions.jl"
        # implement this automatically, generating an array in which the gates are well-ordered
        secondorder_STgates = gen_gates(N, sites, dtau = dtau)
        psi = ST_step(psi, secondorder_STgates)


        if mod(step, 10) == 0
            
            # Measure density profile
            density = [Float64 for _ in 1:N]
            for pos in 1:N
                site_density = op("Nup", sites[pos])
                density[pos] = scalar(imps_expect(psi, site_density, pos))
            end

            plot(density, ylims = (0,1), legend=false,title="Iteració = $step")
            frame(anim)

        end


    end



    gate1 = op("Ntot", sites[j])
    expect2 = scalar(imps_expect(psi, gate1, j))
    #expect2 = imps_expect(psi, gate1, 2)
    print("Escalar final = ", expect2, "\n")
    #print("Escalar final = ", inds(expect2), "\n")
    print()

    # Observable measurement
    #n = expect(psi, "Nup")

    # Generate density plot
    #plot(n, ylimits = (0, 1))
    #savefig("density_1part.png")

    # Generate iTEBD animation
    gif(anim, "1part_evol.gif", fps=6)

    # Generate energy evolution plot
    #plot(energies, yformatter = :scientific, ylimits=(-1e-15, 1e-15))
    #savefig("energy_evol_itebd.png")

end