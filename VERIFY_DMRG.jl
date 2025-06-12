using ITensors, ITensorMPS
#using Revise

# Define a Hamiltonian (replace with your actual Hamiltonian)
# Example: Heisenberg model
N = 10 # Number of sites
t = 1.0
U = 0.0
V = 0.0


function make_hamiltonian(sites; t=1.0, U=8.0, V = 0.0)

    N = length(sites)

    ampo = AutoMPO()
    for i in 1:(N-1)
        # Tunnelling terms
        ampo += -t, "Cdagup", i, "Cup", i+1
        ampo += -t, "Cdagup", i+1, "Cup", i
        ampo += -t, "Cdagdn", i, "Cdn", i+1
        ampo += -t, "Cdagdn", i+1, "Cdn", i
        
        # NN interaction
        ampo += V, "Ntot", i, "Ntot", i+1

        # On-site interaction
        ampo += U, "Nupdn", i
    end

    ampo += U, "Nupdn", N

    H = MPO(ampo, sites)

    return H
end


sites = siteinds("Electron", N; conserve_qns=true)
# Initialize the MPS in a product state
state = [isodd(i) ? "Up" : "Dn" for i in 1:N]
#state[1] = "UpDn"
psi0 = productMPS(sites, state)

H = make_hamiltonian(sites, t=t, U=U, V=V)

# Create an initial random MPS
nsweeps = 5
maxdim = [10, 20, 100, 100, 200]
cutoff = 1e-8

# Run DMRG to find the lowest energy and corresponding MPS
energy,psi = dmrg(H,psi0;nsweeps,maxdim,cutoff)


# Print energy
println("Final energy: ", energy/N)


#error("Para aqui")

# =========================



function make_hamiltonian_gate(sites, j, dt, t, U)
    N = length(sites)
    
    s1 = sites[j]
    s2 = sites[j+1]
    s1p = prime(copy(s1))
    s2p = prime(copy(s2))

    # Kinetic terms
    Hj =  ITensor(dag(s1), dag(s2), s1p, s2p)
    #println("AQUIIII: ", inds(Hj))
    for spin in ("up", "dn")
        c1dag = op("Cdag$spin", s1)
        c1    = op("C$spin", s1)
        c2dag = op("Cdag$spin", s2)
        c2    = op("C$spin", s2)

        Hj += -t * c1dag * c2
        Hj += -t * c2dag * c1
    end

    # On-site interaction terms
    nupdn1 = op("Nupdn", s1) * op("Id", s2)
    nupdn2 = op("Id", s1) * op("Nupdn", s2)
    Hj += U * nupdn1

    if j == (N-1)
        Hj += U * nupdn2
    end

    # Time evolution operator
    expHj = exp(- dt * Hj)
    return Hj, expHj
end




function MPS_energy(psi, H_ops)

    N = length(psi)
    energy = 0

    for j in 1:(N-1)
        psidag_j = dag(prime(psi[j], "Site")) * dag(prime(psi[j+1], "Site"))
        energy += scalar(psidag_j * H_ops[j] * psi[j] * psi[j+1])
    end

    return energy
end

dt = 0.01





# Returns vector of spin correlations, to observe Mott insulating phase
function spin_correlator(psi::MPS, i::Int64, sites::Vector{Index{Vector{Pair{QN, Int64}}}})

    N = length(sites)
    correlations = zeros(N)

    # DEFINE CORRELATOR AS MPO
    
    for j in 1:N
        nextpos = mod1(j + i, N)
        
        ampo = AutoMPO()
        ampo += 1, "Sz", i, "Sz", nextpos
        ampo += 0.5, "S+", i, "S-", nextpos
        ampo += 0.5, "S-", i, "S+", nextpos

        correlator = MPO(ampo, sites)

        corr_spin = inner(psi, apply(correlator, psi))

        correlations[nextpos] = corr_spin
    end

    return correlations
end



#provasz = spin_correlator(psi, 1, sites)
#println("SiSj DMRG: ", provasz)
#provasz = expect(psi, "Nupdn")
#println("NUPDN DMRG: ", provasz)