using ITensors, ITensorMPS




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
    lower_tftensor = dag(copy(transfer_tensor))

    ## Lateral virtual index identification
    lind = inds(unitcell[1])[1]
    rind = inds(unitcell[N])[3]

    prime!(lower_tftensor, lind, rind)
    # println("Mirror tf tensor indices: ", inds(lower_tftensor))

    # Contraction over all physical indices at once
    transfer_tensor *= lower_tftensor
    #println("Mirem aqui: ", inds(transfer_tensor))

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

    counter = 0

    while (norm(lvec - lvec_prev) > tolerance)
        lvec_prev = lvec
        lvec *= transfermatrix
        lvec /= norm(lvec)

        lvec = replaceinds(lvec, (other => connect))

        counter += 1
        if (mod(counter, 10000) == 0)
            lvec = randomITensor(connect)
            println("\t > Reset lateral env vector")
            counter = 0
        end
    end

    #println("Computed lateral vector in ", counter, " iterations")
    #separator = combiner(lateralind, prime(lateralind))
    #separator = replaceinds(separator, (inds(separator)[1] => connect))

    #lenv = lvec * separator

    return lvec # lenv
end





N = 2  # Unit cell size
sites = siteinds("Electron", N) # Using Qubit sites for simplicity
bdim = 5 # Bond dimension

# 2. Initialize the iMPS as a Vector{ITensor} with random tensors
links = [Index(bdim, "Link,n=$i") for i in 0:N]
A = [randomITensor(links[i], sites[i], links[i+1]) for i in 1:N]
A = normalise_unitcell(A)

# Find transfer matrix
tfmatrix = transfermatrix(A)
tolerance = 1e-8
lind = inds(A[1])[1]
#println("Tf matrix indices: ", inds(tfmatrix))


# First attempt: left eigenvector found with power method
@time lvec1 = lateral_env(tfmatrix, tolerance, lind)

tfmatrix2 = matrix(tfmatrix)
#println("Indices tfmatrix: ", inds(tfmatrix))
# Second attempt:
@time eigen_results = eigen(tfmatrix2)
eigenvals = eigen_results.values
eigenvecs = real(eigen_results.vectors)
maxeigval_index = argmax(abs.(eigenvals))
lvec2 = eigenvecs[:,maxeigval_index]
#println("AAAA", lvec2)

# Trobem la diferencia
lv1 = inds(lvec1)
norma = 0.0

for k in 1:length(lvec2)
#    println(lvec1[k])
    norma += abs(lvec1[k] - lvec2[k])
end


itf1, itf2 = inds(tfmatrix)
provemho1 = lvec1 * tfmatrix
provemho1 = replaceinds(provemho1, (itf2 => itf1))
provemho2 = tfmatrix2 * lvec2


println("Diferencia: ", norma)
println("Diferencia primer metode: ", norm(provemho1 - lvec1))
println("Diferencia segon metode: ", norm(provemho2 - lvec2))
println("Mod of lvec1: ", norm(lvec1))
println("Mod of lvec2: ", norm(lvec2))