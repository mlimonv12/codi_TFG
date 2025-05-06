using ITensors

# Define physical index ("site") space
d = 4 # e.g., "Electron" site: empty, up, down, up+down
sA = Index(d, "sA")
sB = Index(d, "sB")
l = Index(10, "link") # initial bond dimension

# Initial random tensors A, B
A = randomITensor(l', sA, l)
B = randomITensor(l', sB, l)

# Define Hamiltonian two-site gate
function make_gate(s1, s2; dt=0.01)
    # example: simple hopping for now
    H = OpSum()
    H += "Cdagup", 1, "Cup", 2
    H += "Cdagup", 2, "Cup", 1
    #H2 = MPO(H, siteinds = [s1, s2])
    # for Hubbard, you'd add on-site terms here
    return exp(-1im * dt * H)
end

gate = make_gate(sA, sB)

# iTEBD sweep parameters
steps = 1000
cutoff = 1e-8
maxdim = 100

for step in 1:steps
    # 1. Form two-site tensor
    theta = A * B

    # 2. Apply gate
    theta = gate * theta

    # 3. SVD and truncate
    U,S,V = svd(theta, (sA,))
    newl = commonind(U, S)
    truncate!(S; cutoff=cutoff, maxdim=maxdim)

    # 4. Normalize
    normS = norm(S)
    S ./= normS

    # 5. Update tensors
    A = U * S
    B = V
end
