using ITensors, ITensorMPS
using Plots
gr()


anim = Animation()

# Codi bàsic per estudiar com de correcte és iTEBD. Generem una partícula tancada
# a una caixa amb μ=0, n=1. Esperem obtenir una funció sinus

# Generates a FH Hamiltonian MPO for U, μ values
function FH_Hamiltonian(N, sites, U, μ)
    t = 1.0 # Hopping amplitude
    os = OpSum()
    for j in 1:N

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


let 
    # Create 10 dim=4 indices (electrons, fermions with S=1/2)
    N = 5
    sites = siteinds("Electron", N)

    # Define simulation parameters
    U = 0.0
    μ = 0.0

    Htotal = FH_Hamiltonian(N, sites, U, μ)

    # Generate MPS with one single particle
    state = fill("Emp", N)
    state[1] = "Up" # Introduce particle
    psi = productMPS(sites, state)
    print(psi)
    psi0 = psi
    #psi = randomMPS(sites,N)
    
    # iTEBD parameters, gates and operator
    cutoff = 1e-5
    maxdim = 64
    dtau = 0.01
    steps = 1500

    # bonds = [(1,2)]  # NN bonds
    # gates = []

    # Suzuki-Trotter decomposition (even-odd gate ordering)
    # Reduces TEBD error


    even_gates = []
    odd_gates = []

    # for (i,j) in bonds
    # Generate hamiltonian gate array
    for i in 1:N

        # Periodic boundary conditions
        if i==N
            iplusone = 1
        else
            iplusone = i + 1
        end

        # Operator definition
        Cup1 = op("Cup", sites[i])
        Cdagup1 = op("Cdagup", sites[i])
        Cdn1 = op("Cdn", sites[i])
        Cdagdn1 = op("Cdagdn", sites[i])
        Cup2 = op("Cup", sites[iplusone])
        Cdagup2 = op("Cdagup", sites[iplusone])
        Cdn2 = op("Cdn", sites[iplusone])
        Cdagdn2 = op("Cdagdn", sites[iplusone])

        # Build H manually
        H = Cup1 * Cdagup2 + Cdagup1 * Cup2# + Cdn1 * Cdagdn2 + Cdagdn1 * Cdn2

        gate = exp(-dtau * H)

        # Since bond pairs are ordered in ascending magnitude, min(i,j) can serve
        # to identify bonds in odd/even positions
        if isodd(i)
            push!(odd_gates, gate)
        else
            push!(even_gates, gate)
        end

    end
    
    #odd_gates, even_gates = get_evenodd_gates(sites, bonds)

    energies = []

    # Run iTEBD algorithm
    for step in 1:steps

        for gate in odd_gates # Apply each gate individually
            psi = apply(gate, psi; cutoff = cutoff)
        end

        for gate in even_gates # Apply each gate individually
            psi = apply(gate, psi; cutoff = cutoff)
        end

        normalize!(psi) # Remember that iTEBD is a dissipative method

        if mod(step,100)==0
            energy = inner(psi', Htotal, psi)
            push!(energies, energy)
            print("Iteració: ", step, " \t | Energia = ", energy, "\n")
            
            # Add frame to video
            n = expect(psi, "Nup")
            plot(n, ylims=(0,1),legend=false,title="Iteració = $step")
            frame(anim)
        end
    end

    # Observable measurement
    n = expect(psi, "Nup")

    # Generate density plot
    plot(n, ylimits = (0, 1))
    savefig("density_1part.png")

    # Generate iTEBD animation
    gif(anim, "iTEBD_evol.gif", fps=10)

    # Generate energy evolution plot
    plot(energies, yformatter = :scientific, ylimits=(-1e-15, 1e-15))
    savefig("energy_evol_itebd.png")

end

# Run iTEBD algorithm
for step in 1:steps

    print("Step = ", step, " out of ", steps, "\n")

    # Odd gates, first application
    for gate in odd_gates
        aaa
    end

    #normalize!(psi) # Remember that iTEBD is a dissipative method
    for i in 1:N
        psi[i] /= norm(psi[i])
    end
    # psi /= norm(psi)

    if mod(step,100)==0
        # energy = inner(psi', Htotal, psi)
        # energy = psi * Htotal * psi
        energy = 0
        #for site in 1:N
        #    energy += conj(psi[site]) * total_gates[site] * psi[site]
        #end
        site=1
        #caca1 = conj(psi[site]) * total_gates[site] * psi[site]
        cacaprod = psi[site]*psi[site+1]
        print("indexos cacaprod: ", inds(cacaprod), "\n")
        print("indexos cacaprod: ", inds(cacaprod)[2], "\n")
        print("Indexos gate: ", inds(total_gates[site]), "\n")
        #caca1 = inner(cacaprod, total_gates[site], cacaprod)
        caca1 = prime(cacaprod, [inds(cacaprod)[2], inds(cacaprod)[3]]) * total_gates[site] * cacaprod
        #prime!(cacaprod, [2,3])
        caca2 = total_gates[site]
        print("caca escalar prova: ", caca1[1], "\n")
        

        push!(energies, energy)
        print("Iteració: ", step, " \t | Energia = ", energy, "\n")
        
        # Add frame to video
        #n = expect(psi, "Nup")
        #plot(n, ylims=(0,1),legend=false,title="Iteració = $step")
        #frame(anim)
    end
end