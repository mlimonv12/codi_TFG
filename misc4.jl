using ITensors, ITensorMPS

# Define MPS
N = 2 # Unit cell size
bdim = 5
sites = siteinds("Electron", N)
#sites = [Index(4, "link-$i") for i in 1:N]
links = [Index(bdim, "link-$i") for i in 0:N]

# Generate iMPS explicitly
psi = [randomITensor(links[i], sites[i], links[i+1]) for i in 1:N]

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
    normal_unitcell = [site/sqrt(dominant_eigenvalue) for site in unitcell]

    return normal_unitcell, dominant_eigenvalue
end


# TO COMPUTE EXPECTATION VALUES OF THE MPS:

# 1. LEFT-CANONICALISE UNIT CELL
# 2. COMPUTE LEFT ENV.
# 3. RIGHT-CANONICALISE UNIT CELL (ORIGINAL COPY)
# 4. COMPUTE RIGHT ENV.
# APPLY BOTH ENVS.






function left_env(tfmatrix::ITensor)

    # Diagonalise transfer matrix


    #


    return lenv
end



println("Comença")

prova1 = transfermatrix(psi)
println(typeof(prova1))
prova2 = normalise_unitcell(psi)

println("Acaba")
println(inds(prova1))