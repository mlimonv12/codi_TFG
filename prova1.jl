using ITensors, ITensorMPS
using Plots
gr()

#include("functions.jl") # Removed this line as it's not needed for the example to run

# Generates a FH Hamiltonian MPO for U, μ values
function FH_Hamiltonian(N, sites, U, μ)
    t = 1.0 # Hopping amplitude
    os = OpSum()
    for j in 1:(N - 1)

        # Hopping terms
        os += -t, "Cdagup", j, "Cup", j + 1 # Changed mod1(j + 1, N) to j+1
        os += -t, "Cup", j, "Cdagup", j + 1 # Changed mod1(j + 1, N) to j+1
        os += -t, "Cdagdn", j, "Cdn", j + 1 # Changed mod1(j + 1, N) to j+1
        os += -t, "Cdn", j, "Cdagdn", j + 1 # Changed mod1(j + 1, N) to j+1

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

    psi = [randomITensor(links[i], sites[i], links[i+1]) for i in 1:N]

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
        T = ITensor(links[i], sites[i], links[i+1])
        T[links[i]=>bdim, sites[i]=>references[state[i]], links[i+1]=>bdim] = 1.0
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
function gen_gates(N::Int64, sites::Vector{Index{Int64}}, dtau::Float64; secondorder = true)

    gates = [ITensor() for _ in 1:N] # Define gates as an empty Vector{ITensor}
    
    for i in 1:N

        # Periodic boundary conditions.  The original code had mod1 here, which is correct
        # for a *periodic* system.  However, the main part of the code was setting up an
        # *open* system.  I've changed this to i+1, and added a check to prevent indexing
        # beyond the end of the sites array.
        iplusone = (i < N) ? i + 1 : 1

        # h = get_operator("TunnellingUP", sites, i)
        h = op("Cup", sites[i]) * op("Cdagup", sites[iplusone])
        h += op("Cdagup", sites[i]) * op("Cup", sites[iplusone])
        h *= -1

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
        #println(inds(gate))

        gates[i] = gate 
    end

    return gates
end


# Applies a two-site gate to an iMPS of type Vector{ITensor}, and returns the split iMPS sites
function apply_gate(site1::ITensor, site2::ITensor, gate::ITensor; cutoff = 1e-20, maxdim = 30)

    # Check for double index connection: this is the case of a looped two-site unit cell iMPS
    l1, p1, r1 = inds(site1)
    l2, p2, r2 = inds(site2)

    #println("SITE1: ", inds(site1))
    #println("SITE2: ", inds(site2))

    # Prime leftmost internal index to prevent double contractions
    prime!(site1, l1)

    # If the double connection is given, one of the lateral indices is primed to avoid double contractions
    #if (l1 == r2) && (r1 == l2)
    #   println("caca aqui")
    #   prime!(site1, l1)
    #end

    # Contract the site at pos with the next one, two-site gate
    prodsite = site1 * site2

    # Contract the joint site with the two-site operator
    prodsite_evol = gate * prodsite
    #println("123", inds(prodsite_evol))

    # Define combiner tensors for left and right sides: this will help switch prodsite_evol
    # from a D,d,d,D indexed tensor to a D*d, D*d tensor, enabling SVD
    i1, i2, i3, i4 = inds(prodsite_evol)
    cleft = combiner(i1, i3)
    cright = combiner(i2, i4)

    # Apply combiners, making prodsite_matrix a D*d, D*d tensor
    prodsite_matrix = prodsite_evol * cleft
    prodsite_matrix *= cright

    # Perform SVD
    U, S, V = svd(prodsite_matrix, inds(prodsite_matrix)[1]; cutoff = cutoff, maxdim = maxdim) # Added cutoff

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
    replaceindex!(sql, inds(sql)[2], r1)
    replaceindex!(sqr, inds(sqr)[1], r1)
    
    # Finally form and reshape the new site tensors, of indices D,d,D
    #println("inds U: ", inds(U))
    #println("inds sql: ", inds(sql))
    #println("inds cleft: ", inds(cleft))
    L = U * sql * dag(cleft)
    R = sqr * V * dag(cright)
    
    L = permute(L, p1, l1, r1) # Changed the order of indices to match
    L = noprime(L)
    R = noprime(R)
    #println("Indexos de L: ", inds(L))
    #println("Indexos de R: ", inds(R), "\n")

    # Return updated sites
    return L, R
end


# Applies a Suzuki-Trotter step to an MPS
function ST_step(iMPS::Vector{ITensor}, gates::Vector{ITensor}; secondorder = true, maxdim = 30, cutoff = 1e-20) # Added cutoff

    N = length(iMPS)

    # Odd sites
    for i in 1:N
        if isodd(i)
            iMPS[i], iMPS[mod1(i+1, N)] = apply_gate(iMPS[i], iMPS[mod1(i+1, N)], gates[i], maxdim = maxdim, cutoff = cutoff) # Added cutoff
        end
    end
    
    # Even sites
    for i in 1:N
        if iseven(i)
            iMPS[i], iMPS[mod1(i+1, N)] = apply_gate(iMPS[i], iMPS[mod1(i+1, N)], gates[i], maxdim = maxdim, cutoff = cutoff) # Added cutoff
        end
    end
    
    # In second-order ST evolution odd gates are applied again
    if secondorder
        for i in 1:N
            if isodd(i)
                iMPS[i], iMPS[mod1(i+1, N)] = apply_gate(iMPS[i], iMPS[mod1(i+1, N)], gates[i], maxdim = maxdim, cutoff = cutoff) # Added cutoff
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
        operator = op("Cdagup", sites[site]) * op("Cup", sites[nsite])
        operator += op("Cup", sites[site]) * op("Cdagup", sites[nsite])
        operator *= -t

    elseif name == "TunnellingDOWN"
        operator = op("Cdagdn", sites[site]) * op("Cdn", sites[nsite])
        operator += op("Cdn", sites[site]) * op("Cdagdn", sites[nsite])
        operator *= -t

    elseif name == "Tunnelling"
        operator = op("Cdagup", sites[site]) * op("Cup", sites[nsite])
        operator += op("Cup", sites[site]) * op("Cdagup", sites[nsite])
        operator += op("Cdagdn", sites[site]) * op("Cdn", sites[nsite])
        operator += op("Cdn", sites[site]) * op("Cdagdn", sites[nsite])
        operator *= -t
    
    elseif name == "Onsite"
        operator = -U * op("Nupdn", sites[site])

    elseif name == "ChemPot"
        operator = -μ * op("Ntot", sites[site])
    
    elseif name == "Hubbard"
        operator = op("Cdagup", sites[site])*op("Cup", sites[nsite])
        operator += op("Cup", sites[site])*op("Cdagup", sites[nsite])
        operator += op("Cdagdn", sites[site])*op("Cdn", sites[nsite])
        operator += op("Cdn", sites[site])*op("Cdagdn", sites[nsite])
        operator *= -t
        operator += -U * op("Nupdn", sites[site])
    
    else
        error("\n\n No operator identified with the following name: ", name)
        #throw(InterruptException())
    end

    return operator
end


# Returns the expected value of a given operator for an iMPS made of Vector{ITensor}
# Works for operators of any size
# Computes the expected value of an operator in an iMPS Vector{ITensor} object
function imps_expect2(iMPS::Vector{ITensor}, operator::ITensor, site::Int64)

    # Find number of sites the operator acts on
    s = Int64(length(size(operator))/2)
    #print("\n\t> S = ", s, "\n")

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
    applied_tosite = site_tensor' * operator * site_tensor # Corrected order of multiplication

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

    return scalar(expect) # Return a scalar
end


# Canonicalises and normalises the iMPS


# DONE
function normalise_unitcell(unitcell::Vector{ITensor})

    # Get transfer matrix
    tf_matrix = transfermatrix(unitcell)
    linds, rinds = inds(tf_matrix)


    # Diagonalise
    eigenvals, eigenvecs = eigen(tf_matrix, linds, rinds)

    # Get dominant eigenvalue
    dominant_eigenvalue = maximum(abs, eigenvals)

    # Normalise
    N = length(unitcell)
    normal_unitcell = unitcell / sqrt(dominant_eigenvalue^(1/N))

    return normal_unitcell#, dominant_eigenvalue
end



# DONE?
function imps_expect(unitcell::Vector{ITensor}, operator::ITensor; tolerance = 1e-8)

    N = length(unitcell)

    # Define lateralmost indices
    lind = inds(unitcell[1])[1]
    rind = inds(unitcell[N])[3]

    tfmatrix = transfermatrix(unitcell)
    # println("Tf matrix type: ", typeof(tfmatrix), inds(tfmatrix))

    lenv = lateral_env(tfmatrix, tolerance, lind)
    renv = lateral_env(tfmatrix, tolerance, rind; left = false)

    # Contract iMPS cells
    A = unitcell[1]
    for i in 2:N
        A *= unitcell[i]
    end

    Adag = prime(dag(A))

    # println("A initial indices: ", inds(A))
    # println("Operator indices: ", inds(operator))

    # Compute A adjoint
    # Aadj = 

    # Contract R into A
    A = A * renv # Changed order
    # println("AR indices: ", inds(A))

    # Contract At into RA
    A = Adag * A # Changed order
    # println("AtRA indices: ", inds(A))

    # Contract operator into AtRA
    A = operator * A # Changed order
    #println("AtRA + operator indices: ", inds(A))

    # Contract L into final tensor, scalar result
    A = lenv * A # Changed order

    #error("Parem per ara")
    # Return observable
    return scalar(A)
end


# DONE
function transfermatrix(unitcell::Vector{ITensor})

    N = length(unitcell)
    
    # Contract over internal indices
    transfer_tensor = unitcell[1]
    for i in 2:N
        transfer_tensor *= unitcell[i]
    end

    #println("First tf tensor indices: ", inds(transfer_tensor))

    # Contract over physical indices. We prime lateral indices to avoid accidental contraction
    lower_tftensor = copy(transfer_tensor)

    ## Lateral virtual index identification
    lind = inds(unitcell[1])[1]
    rind = inds(unitcell[N])[3]

    prime!(lower_tftensor, lind, rind)
    # println("Mirror tf tensor indices: ", inds(lower_tftensor))

    # Contraction over all physical indices at once
    transfer_tensor *= dag(prime(lower_tftensor, siteinds(unitcell))) # Corrected this line

    # Reshape transfer tensor to matrix
    # println("Transfer tensor indices: ", inds(transfer_tensor), "\n")
    i1, i2, i3, i4 = inds(transfer_tensor)
    cleft = combiner(i1, i3)
    cright = combiner(i2, i4)

    return transfer_tensor * cleft * cright
end


# FIRST ATTEMPT: iterative method
function lateral_env(transfermatrix::ITensor, tolerance::Float64, lateralind::Index{Int64}; left = true)

    l, r = inds(transfermatrix)

    if left
        connect = l
        other = r
    else
        connect = r
        other = l
    end

    lvec = randomITensor(connect)
    lvec_prev = randomITensor(connect)

    while (norm(lvec - lvec_prev) > tolerance)
        lvec_prev = lvec
        lvec = transfermatrix * lvec # Changed order

        lvec /= norm(lvec)

        lvec = replaceinds(lvec, (other => connect))
    end

    separator = combiner(lateralind, prime(lateralind))
    separator = replaceinds(separator, (inds(separator)[1] => connect))

    lenv = lvec * separator

    return lenv
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





# Example usage:
# Assume you have a Vector{ITensor} called 'my_sites'
# my_normalized_sites = normalize_mps!(deepcopy(my_sites))


# Computes the Von Neumann entropy of an iMPS at a given site
function site_entropy(iMPS::Vector{ITensor}, site::Int64)

    

    entropy = 0


    return entropy
end




let 
    # 
    N = 5 # Changed to 5 for this example
    bdim = 16
    sites = siteinds("Electron", N)
    links = [Index(bdim, "link-$i") for i in 0:N]

    # Create an iMPS with a known initial state
    psi = init_iMPS(N, sites, links)
    #print(typeof(sites), "\n")
    psi = set_FHstate(psi, sites, links, [0,1,0,0,0]) # Corrected state vector length

    # Define simulation parameters
    U = 0.0
    μ = 0.0
    
    # iTEBD parameters, gates and operator
    cutoff = 1e-5
    dtau = 0.01
    steps = 100
    checkevery = 10
    framespersecond = 6

    mesures = []
    #println("PRINCIPI ", inds(psi[1]))

    # Generate iTEBD gates
    secondorder_STgates = gen_gates(N, sites, dtau)

    # Prepare iMPS to apply iTEBD
    psi = loop_iMPS(psi)

    anim = Animation() # Define anim here

    # Begin iTEBD loop
    for step in 1:steps

        println("Step ", step, " of ", steps)
        
        # Evolve the system using iTEBD, second-order Suzuki-Trotter gate ordering:
        # First (square-rooted) odd gates are applied, then even, then odd gates again
        # The functions gen_gates and secondorder_STstep defined in the file "functions.jl"
        # implement this automatically, generating an array in which the gates are well-ordered
        psi = ST_step(psi, secondorder_STgates, maxdim = bdim, cutoff = cutoff) # Added cutoff

        # Un-loop iMPS
        psi = unloop_iMPS(psi, links[1])

        # Normalise iMPS
        psi = normalise_unitcell(psi)

        # Re-loop iMPS
        psi = loop_iMPS(psi)


        if mod(step, checkevery) == 0

            # Unloop iMPS to compute expectation values
            psi = unloop_iMPS(psi, links[1])
            
            # Measure density profile
            density = zeros(N)
            for pos in 1:N
                site_density = op("Nup", sites[pos])
                density[pos] = imps_expect(psi, site_density) # Corrected the call to imps_expect
            end

            plot(density, ylims = (0,1), legend=false, title="Iteration = $step")
            frame(anim)

            # Re-loop iMPS to continue iTEBD
            psi = loop_iMPS(psi)

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
    gif(anim, "1part_evol.gif", fps=framespersecond)

    # Generate energy evolution plot
    #plot(energies, yformatter = :scientific, ylimits=(-1e-15, 1e-15))
    #savefig("energy_evol_itebd.png")

end
