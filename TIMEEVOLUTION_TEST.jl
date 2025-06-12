using ITensors, ITensorMPS
using Plots, Statistics
using HDF5, DelimitedFiles   # Data exportation
gr()

anim = Animation()



function make_hamiltonian(sites, t, U; dt=0)

    N = length(sites)

    ampo = AutoMPO()
    for i in 1:(N-1)
        # Tunnelling terms
        ampo += -t, "Cdagup", i, "Cup", i+1
        ampo += -t, "Cdagup", i+1, "Cup", i
        ampo += -t, "Cdagdn", i, "Cdn", i+1
        ampo += -t, "Cdagdn", i+1, "Cdn", i

        # On-site interaction
        ampo += U, "Nupdn", i
    end

    ampo += U, "Nupdn", N

    H = MPO(ampo, sites)

    item = H[1]

    for i in 2:N
        item *= H[i]
    end

    expH = exp(-dt * item)

    return H, expH
end


function make_hamiltonian_gate(sites, j, dt, t, U; secondorder = false)
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

    if (isodd(j) && secondorder)
        expHj = exp(- dt/2 * Hj)
    else
        expHj = exp(- dt * Hj)
    end

    return Hj, expHj
end



function make_hamiltonian_gate_inspired(H, sites, j, dt, t, U)

    Hj = H[j] * H[j+1]

    return Hj
end



function make_tunnelling_gate(sites, j, dt, t)
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

    # Time evolution operator
    expHj = exp(- dt * Hj)

    return expHj
end


function make_onsite_gate(sites, j, dt, U)
    N = length(sites)
    s1 = sites[j]
    s1p = prime(copy(s1))

    # Kinetic terms
    Hj =  ITensor(dag(s1), s1p)
    
    # On-site interaction terms
    Hj += U * op("Nupdn", s1)

    # Time evolution operator
    expHj = exp(- dt * Hj)

    return expHj
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



# Returns vector of spin correlations, to observe Mott insulating phase
function spin_correlator(psi::MPS, N::Int64, site::Int64, sites)

    correlations = zeros(N)
    
    for j in 1:(N-1)
        nextpos = mod1(j + site, N)
        psidag_ij = dag(prime(psi[site], "Site")) * dag(prime(psi[nextpos], "Site"))
        SiSj = op(sites, "Sz", site) * op(sites, "Sz", nextpos)
        corr_spin = scalar(psidag_ij * SiSj * psi[site] * psi[nextpos]) * 2
        
        correlations[j] = corr_spin
    end

    return correlations
end







# Define the number of sites
N = 4

if mod(N,2) != 0
    error("N ha de ser parell")
end



# Create site indices for the system
sites = siteinds("Electron", N)#; conserve_qns=true)

# Initialize the MPS in a product state
state = [isodd(i) ? "Up" : "Dn" for i in 1:N]
psi0 = productMPS(sites, state)
psi = copy(psi0)

global wavefunction = psi[1]
for j in 2:N
    wavefunction *= psi[j]
end



# Define the Hamiltonian parameters
t = 1.0
U = 4.0
savedata = false
steps = 3000#600
checkevery = 50#50
framespersecond = 10
cutoff = 1e-8
bdim = 20
# Define the time step
dt = 0.01
saveconvergence = true



println("Inicia programa")



H, expH = make_hamiltonian(sites, t, U, dt = dt)

Hmatrix = H[1]
for j in 2:N
    Hmatrix *= H[j]
end


# Make hamiltonian & TEBD gates
H_ops = []
gates = []
tunnelgates = []
onsitegates = []
energia = 0.0



densupdn = []
energies = [[], []]

println("Comencem a aplicar gates")
# Apply gates
for step in 1:steps
    #end

    wavefunction *= expH
    normalize!(wavefunction)

    if mod(step,checkevery)==0
        println("Step ", step, " of ", steps)
        densup = expect(psi, "Nup")
        densdn = expect(psi, "Ndn")
        densupdn = expect(psi, "Nupdn")
        plot(densup, ylims = (0,1), lc=:red)
        plot!(densdn, ylims = (0,1), lc=:blue)
        plot!(densupdn, ylims = (0,1), lc=:green)
        frame(anim)
    
        
        #energia = inner(psi, apply(H, psi))
        #energia = inner(dag(psi), H, psi)
        #global energia = MPS_energy(psi, H_ops)
        energy = scalar(wavefunction * noprime(Hmatrix * wavefunction))

        println("ENERGIA MANUAL: ", energy)
        push!(energies[1], step)
        push!(energies[2], energy)

        #push!(energies, energia)
    end
end

if savedata
    #prelude = "figures/periodic/op_"*operator*"/"
    #fig_id = "U_"*string(U)*"_V_"*string(V)*"_dt_"*string(dtau)*"_bdim_"*string(bdim)*"_steps"*string(steps)
    gif(anim, "provaopen1.gif", fps=framespersecond)

    plot(energies)
    savefig("Energia_evol_OBC.png")
end



#println("Valor final: ", mean(densupdn), " ± ", std(densupdn))
println("Energia final: ", energy)

println("Fi del programa")


title = "N_$N-U_$U-V_$V-bdim_$bdim-dt_$dt-o2_zero"

# Save to file
if saveconvergence
    writedlm( "converged/Energy/$title.csv", energies, ',')
end
