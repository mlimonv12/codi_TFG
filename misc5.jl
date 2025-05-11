using ITensors, ITensorMPS

N = 2  # Unit cell size
sites = siteinds("Qubit", N) # Using Qubit sites for simplicity
d = 2    # Dimension of the local Hilbert space
bdim = 4 # Bond dimension

# 2. Initialize the iMPS as a Vector{ITensor} with random tensors
links = [Index(bdim, "Link,n=$i") for i in 0:N]
A = [randomITensor(links[i], sites[i], links[i+1]) for i in 1:N]

# For an infinite MPS, the first and last link index are the same
A[1] = randomITensor(links[1], sites[1], links[2])
A[2] = randomITensor(links[2], sites[2], links[3])

# 3. Canonicalize the iMPS (Mixed Canonical Form)
#   ITensors.jl doesn't have a direct function to canonicalize a Vector{ITensor}
#   representing an iMPS.  Canonicalization is more naturally done with the InfiniteMPS type.
#   However, we can bring it to a kind of left/right canonical form.

function canonicalize_iMPS!(A::Vector{ITensor})
    N = length(A)
    for i in 1:N
        # Bring A[i] to right canonical form
        if i < N
            Q, R = qr(A[i], [inds(A[i])[1], inds(A[i])[2]]); # Corrected qr
            A[i] = Q
            A[i+1] = R * A[i+1]
        else
            Q, R = qr(A[i], [inds(A[i])[1], inds(A[i])[2]]); # Corrected qr
            A[i] = Q
            #A[1] = R * A[1] # close the loop.  This is necessary for a true iMPS.
        end
    end
    for i in N:-1:1
            if i > 1
            L,P = lq(A[i], [inds(A[i])[2], inds(A[i])[3]]); # Corrected lq
            A[i] = P
            A[i-1] = A[i-1] * L
            else
                L,P = lq(A[i], [inds(A[i])[2], inds(A[i])[3]]); # Corrected lq
                A[i] = P
                #A[N] = A[N] * L
            end
    end
end

canonicalize_iMPS!(A)

println("Inds A:", inds(A))