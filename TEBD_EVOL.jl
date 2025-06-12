using ITensors, ITensorMPS  # Tensor elements
using Plots, Statistics     # Data processing
using HDF5, DelimitedFiles  # Data exportation
gr()

#@set_warn_order 30

anim = Animation()



# Makes hamiltonian MPO
function make_hamiltonian(sites, t, U, V)

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

    #item = H[1]

    #for i in 2:N
    #    item *= H[i]
    #end

    return H
end


# Makes gates for iTEBD
function make_hamiltonian_gate(sites, j, t, U, V, dt; secondorder = false)

    N  = length(sites)

    ampo = AutoMPO()
    # Tunnelling terms
    ampo += -t, "Cdagup", 1, "Cup", 2
    ampo += -t, "Cdagup", 2, "Cup", 1
    ampo += -t, "Cdagdn", 1, "Cdn", 2
    ampo += -t, "Cdagdn", 2, "Cdn", 1
    
    # NN interaction
    ampo += V, "Ntot", 1, "Ntot", 2

    # On-site interaction
    ampo += U, "Nupdn", 1

    if (j == (N-1))
        ampo += U, "Nupdn", 2
    end
        
    Hj = MPO(ampo, sites[j:j+1])
    Hj_comb = Hj[1] * Hj[2]

    # Time evolution operator
    if (isodd(j) && secondorder)
        expHj = exp(- dt/2 * Hj_comb)
    else
        expHj = exp(- dt * Hj_comb)
    end

    # Fill hamiltonian gate with identities to allow for proper MPO application
    #for k in 3:N
    #    nextpos = mod1(j + k, N)
    #    idsite = AutoMPO()
    #    idsite += 1, "Id", 1
    #    push!(Hj, MPO(idsite, [sites[nextpos]]))
    #end

    return Hj, expHj
end



function MPS_energy(psi, H_ops)

    N = length(psi)
    energy = 0

    for j in 1:(N-1)
        contracted = H_ops[j][1] * H_ops[j][2]

        energy += inner(psi, contracted * psi)
        #energy += inner(psi, apply(H_ops[j], psi))
    end

    return energy
end




# Clears a printed line and prints another one
# (source: https://discourse.julialang.org/t/how-clear-the-printed-content-in-terminal-and-print-to-the-same-line/19549/3)
function overprint(str)  
    print("\u1b[1F")
    #Moves cursor to beginning of the line n (default 1) lines up   
    print(str)   #prints the new line
   print("\u1b[0K") 
   # clears  part of the line.
   #If n is 0 (or missing), clear from cursor to the end of the line. 
   #If n is 1, clear from cursor to beginning of the line. 
   #If n is 2, clear entire line. 
   #Cursor position does not change. 

    println() #prints a new line, i really don't like this arcane codes
end



# ====================================================================================================

# CONTROL VARIABLES
function TEBD_EVOL(
    # Number of sites
    N;
    # Hamiltonian parameters
    t = 1.0,
    U = 4.0,
    V = 0.0,
    # Simulation parameters
    dt = 0.01,
    maxsteps = 100,
    checkevery = 20,
    framespersecond = 10,
    cutoff = 1e-8,
    bdim = 20,
    ref_energy = 1000.0,
    energy_precision = 1e-7,
    secondorder = false,
    monitorobservable = true,
    savefigures = true,
    saveconvergence = false,
    saveMPS = true
)

# ====================================================================================================


    println("Simulating N=$N, U=$U, V=$V")


    # Create site indices for the system. QNs are conserved to 
    sites = siteinds("Electron", N; conserve_qns=true)


    # Initialize the MPS in a product state
    state = [isodd(i) ? "Up" : "Dn" for i in 1:N]
    psi0 = productMPS(sites, state)
    psi = copy(psi0)


    # Generate Hamiltonian
    println("Generating Hamiltonian MPO...")
    H = make_hamiltonian(sites, t, U, V)
    overprint("Generated hamiltonian MPO")

    # Make TEBD gates
    H_ops = []
    gates = []
    
    println("Generating TEBD gates...")
    for j in 1:(N-1)
        Hj, gate = make_hamiltonian_gate(sites, j, t, U, V, dt, secondorder = secondorder)
        push!(H_ops, Hj)
        push!(gates, gate)
    end
    overprint("Generated TEBD gates")


    # Define variables
    densupdn = []
    energies = [[], []]
    observe = [[], [], []]
    energy = 0.0
    prev_energy = 100.0


    # Begin TEBD
    println("> Start TEBD evolution")
    println()

    for step in 1:maxsteps

        # Odd gates
        for i in 1:(N-1)
            if isodd(i)
                psi = apply(gates[i], psi, cutoff=cutoff, maxdim=bdim)
            end
        end
        normalize!(psi)
        
        # Even gates
        for i in 1:(N-1)
            if iseven(i)
                psi = apply(gates[i], psi, cutoff=cutoff, maxdim=bdim)
            end
        end
        normalize!(psi)

        # Odd gates again, in case of second-order Suzuki-Trotter gate application
        if secondorder
            for i in 1:(N-1)
                if isodd(i)
                    psi = apply(gates[i], psi, cutoff=cutoff, maxdim=bdim)
                end
            end
        end
        normalize!(psi)


        if mod(step,checkevery)==1
            overprint("Step $step of $maxsteps")
        
            # Find energy, ⟨ψ|H|ψ⟩
            #energy = MPS_energy(psi, H_ops)
            energy = inner(psi, apply(H, psi))

            # println("Energy at step $step: ", energy)

            # Check for convergence. If desired we can specify a given reference energy
            if (ref_energy != 1000.0)
                prev_energy = ref_energy
            end

            if (abs(prev_energy - energy) < energy_precision)
                println("Converged energy to a precision of $energy_precision")
                break
            end

            push!(energies[1], step)
            push!(energies[2], energy)


            if monitorobservable
                obsval = expect(psi, "Nupdn")

                push!(observe[1], step)
                push!(observe[2], mean(obsval))
                push!(observe[3], std(obsval))
            end

            if step == maxsteps
                println("Converged energy to a precision of ", abs(energy - prev_energy))
            end

            prev_energy = energy

        end
    end


    if savefigures
        #gif(anim, "provaopen1.gif", fps=framespersecond)

        plot(energies[1], energies[2], lw=4, xlabel = "Iteracions", ylabel = "Energia")
        savefig("Energia_evol_OBC.png")
    end


    println("Final energy: $energy")

    title = "N_$N-U_$U-V_$V-bdim_$bdim-dt_$dt-o2_$secondorder"

    # Save to file
    if saveconvergence
        writedlm( "converged/Energy/$title.csv", energies, ',')

        if monitorobservable
            writedlm( "converged/Other/$title.csv", observe, ',')
        end
    end


    if saveMPS
        f = h5open("converged/MPS/$title.h5","w")
        write(f,"T",psi)
        close(f)
    end


    println("--- End simulation ---")

    return
end

