using ITensors, ITensorMPS

# Define MPS
N = 2 # Unit cell size
bdim = 10
sites = siteinds("Electron", N)
#sites = [Index(4, "link-$i") for i in 1:N]
links = [Index(bdim, "link-$i") for i in 0:N]

# Generate iMPS explicitly
psi = [randomITensor(links[i], sites[i], links[i+1]) for i in 1:N]

# PSEUDOCODE
function expected_imps(unitcell::Vector{ITensor}, operator; tolerance = 1e-8)

    lenv = lateral_env(unitcell)
    renv = lateral_env(unitcell; left = false)

    # Contract iMPS cells
    A = unitcell[1]
    for i in 2:N
        A *= unitcell[i]
    end
    
    # define A'
    lind = inds(unitcell[1])[1]
    rind = inds(unitcell[N])[3]
    #Aprime = prime(A, lind, rind)

    # Contract R into A
    A *= renv

    # Contract A' into RA
    A *= prime(A, lind, rind)

    # Contract operator into A'RA
    A *= operator

    # Contract L into final tensor, scalar result
    A *= lenv

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

    # Contract over physical indices. We prime lateral indices to avoid accidental contraction
    lower_tftensor = copy(transfer_tensor)

    ## Lateral virtual index identification
    lind = inds(unitcell[1])[1]
    rind = inds(unitcell[N])[3]

    prime!(lower_tftensor, lind, rind)

    # Contraction over all physical indices at once
    transfer_tensor *= lower_tftensor

    # Reshape transfer tensor to matrix
    i1, i2, i3, i4 = inds(transfer_tensor)
    cleft = combiner(i1, i3)
    cright = combiner(i2, i4)

    return transfer_tensor * cleft * cright
end

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

    return normal_unitcell, dominant_eigenvalue
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
        lvec *= transfermatrix

        lvec /= norm(lvec)

        lvec = replaceinds(lvec, (other => connect))
    end

    separator = combiner(lateralind, prime(lateralind))
    separator = replaceinds(separator, (inds(separator)[1] => l))

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







println("Comença")


prova2, dom_eig = normalise_unitcell(psi)
println(typeof(prova2))
prova1 = transfermatrix(prova2)

println("Acaba")
println(inds(prova1))

# Check normalisation
lind, rind = inds(prova1)
evals, evecs = eigen(prova1, lind, rind)

#evals_conj = adjoint(evals)

#evals *= evals_conj

eval_vec = [abs(evals[i,i]) for i in 1:size(evals)[1]]

j=1
println("eigenvalues: ", maximum(eval_vec))

