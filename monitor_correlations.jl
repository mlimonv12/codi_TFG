using ITensors, ITensorMPS
using HDF5, Statistics
using Plots, DelimitedFiles
gr()

function get_name(N, U, V, bdim, dt, secondorder)
    name = "converged/MPS/N_$N-U_$U-V_$V-bdim_$bdim-dt_$dt-o2_$secondorder.h5"
    return name
end



# Returns vector of spin correlations, to observe Mott insulating phase
function spin_correlator(psi::MPS, i::Int64)

    sites = siteinds(psi)
    N = length(sites)
    correlations = zeros(N)

    # DEFINE CORRELATOR AS MPO
    println("Enter correlator")
    
    for j in 1:(N)

        k = mod1(i + j, N)
    
        # Define MPO
        ampo = AutoMPO()
        #ampo += 1, "Nup", i, "Nup", k
        #ampo += -1, "Nup", i, "Ndn", k
        #ampo += -1, "Ndn", i, "Nup", k
        #ampo += 1, "Ndn", i, "Ndn", k
        ampo += 1, "Sz", i, "Sz", k
        
        ampo += 0.5, "S+", i, "S-", k
        ampo += 0.5, "S-", i, "S+", k

        sitespin = MPO(ampo, sites)

        #corr_spin = inner(psi, psi * op("Sz", sites[j]))

        corr_spin = inner(psi, apply(sitespin, psi))

        correlations[k] = corr_spin
    end

    return correlations
end


# Returns vector of spin correlations, to observe Mott insulating phase
function spin_correlator_optimal(psi::MPS, i::Int64)

    sites = siteinds(psi)
    N = length(sites)
    correlations = zeros(N)

    # DEFINE CORRELATOR AS MPO
    println("Enter correlator $i")
    
    for j in (i-1):(N)

        k = mod1(j, N)
    
        # Define MPO
        ampo = AutoMPO()
        #ampo += 1, "Nup", i, "Nup", k
        #ampo += -1, "Nup", i, "Ndn", k
        #ampo += -1, "Ndn", i, "Nup", k
        #ampo += 1, "Ndn", i, "Ndn", k
        ampo += 1, "Sz", i, "Sz", k
        
        ampo += 0.5, "S+", i, "S-", k
        ampo += 0.5, "S-", i, "S+", k

        sitespin = MPO(ampo, sites)

        #corr_spin = inner(psi, psi * op("Sz", sites[j]))

        corr_spin = inner(psi, apply(sitespin, psi))

        correlations[k] = corr_spin
    end

    return correlations
end


# Compute the BOW order parameter, to study the effect of increasing V
function BOW_param(psi::MPS)

    sites = siteinds(psi)
    N = length(sites)
    llim = Int64(N/2)
    BOWs = zeros(2)

    # DEFINE AS MPO
    
    for i in llim:(llim+1)
        ampo = AutoMPO()
        ampo += 1, "Cdagup", i, "Cup", i+1
        ampo += 1, "Cdagup", i+1, "Cup", i
        ampo += 1, "Cdagdn", i, "Cdn", i+1
        ampo += 1, "Cdagdn", i+1, "Cdn", i
        #ampo += 1, "Sz", i, "Sz", nextpos
        #ampo += 0.5, "S+", i, "S-", nextpos
        #ampo += 0.5, "S-", i, "S+", nextpos

        BOW_op = MPO(ampo, sites)

        BOW_site = inner(psi, apply(BOW_op, psi))

        BOWs[i-llim+1] = BOW_site
    end

    return abs(BOWs[1] - BOWs[2])
end

# Compute the CDW order parameter for a given MPS
function CDW_param(psi::MPS, N::Int64)

    ens = expect(psi, "Ntot")

    ΔCDW = 0.0
    for i in 1:N
        ΔCDW += (-1)^i*ens[i]
    end

    return ΔCDW
end

# Compute the numerical derivative of a vector
function num_der(xdata, ydata)

    N = length(ydata)
    xnew = []
    ynew = []

    for i in 1:(N-1)
        push!(ynew, (ydata[i+1] - ydata[i])/(xdata[i+1] - xdata[i]))
        push!(xnew, (xdata[i+1] + xdata[i])/2)
    end

    return xnew, ynew
end



N = 10
U = 10.0
V = 0.0
bdim = 20
dt = 0.01
secondorder = false
study_cdw = false



#Vs = [1.9, 1.92, 1.94, 1.96, 1.98, 2.00, 2.02, 2.04, 2.06, 2.08, 2.1, 2.12, 2.14, 2.16]
#Vs = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0]#, 10.0, 12.0, 14.0] # U = 10

#Vs = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0]
Us = [4.0, 6.0, 8.0, 10.0, 15.0, 20.0]
results = [[], [], [], [], [], []]#zeros(length(Vs))
Vs = [[], [], [], [], [], []]

if study_cdw
    #for item in 1:1
    for i in 1:6
        U2 = Us[i]
        Vxs = []

        for j in 0:28
            if j <= 4
                push!(Vxs, round(j*U2*0.1, digits = 1))
            elseif j > 4 && j <= 13
                push!(Vxs, round(((j-4)*0.02 + 0.4)*U2, digits = 2))
            else
                push!(Vxs, round((j-8)*U2*0.1, digits = 1))
            end
        end


    #    if i != 1

        for item in Vxs

            #Vinst = round(Vx*Us[i]*0.1, digits=1)
            #push!(Vxs, Vinst)

            title = get_name(N, Us[i], item, bdim, dt, secondorder)

            # Load an MPS file
            f = h5open(title,"r")
            psi = read(f,"T",MPS)
            close(f)

            cdw1 = CDW_param(psi, N)
            
            push!(results[i], cdw1/U2)
        end

        Vs[i] = Vxs/U2

        println(length(results[i]))

        if i == 1
            plot(Vxs/U2, results[i], xlabel = "V/U", ylabel = "Δ_CDW", label = "U=$U2", xlims = (0,2), lw = 3)
        else
            plot!(Vxs/U2, results[i], xlabel = "V/U", ylabel = "Δ_CDW", label = "U=$U2", xlims = (0,2), lw = 3, legendfontsize = 14)
        end
    end

    writedlm("processed/CDW_vals.csv", results, ',')
    writedlm("processed/CDW_V.csv", Vs, ',')


    savefig("cdw1.png")
    # Return correlator
    #spin_corr = spin_correlator(psi, 1)


    for i in 1:6
        nderx, ndery = num_der(Vs[i], results[i])
        
        U2 = round(Us[i], digits = 0)

        if i == 1
            plot(nderx, ndery/U2, xlabel = "V/U", ylabel = "dΔ/dV", label = "U=$U2", lw = 2)
        else
            plot!(nderx, ndery/U2, xlabel = "V/U", ylabel = "dΔ/dV", label = "U=$U2", lw = 2)
        end
    end

    savefig("cdw_numder.png")
end

# Return BOW
#BOW_psi = BOW_param(psi)

#push!(results, BOW_psi)
#println(mean(BOW_psi))
#end

#spin_corr = abs.(spin_corr)

#println(spin_corr)
#plot(spin_corr, ylim=(-1, 1))


N = 4

title = get_name(N, U, V, bdim, dt, secondorder)

# Load an MPS file
f = h5open(title,"r")
psi = read(f,"T",MPS)
close(f)


spincorrsvec = zeros(N,N)

for a in 1:N
    spin_corr = spin_correlator(psi, a)
    #push!(spincorrsvec, spin_corr)
    for b in 1:N
        spincorrsvec[a,b] = spin_corr[b]
    end
end

writedlm("processed/spincorrs-N_$N-U_$U-V_$V.csv", results, ',')


output_width = 400
output_height = 360


heatmap(spincorrsvec,
    # xlabel = "i",
    # ylabel = "j",
    yflip = true,
    clim = (-1, 1),
    aspect_ratio = :equal,
    xlim = (0.5, 10.5),
    ylim = (0.5, 10.5),
    xticks = (1:10),
    yticks = (1:10),
    size = (output_width, output_height),
    titlefontsize = 14,
    colorbar_title = "⟨S\$_i\$S\$_j\$⟩", # Optional: title for the colorbar
    colormap = :plasma)      # Choose a colormap (e.g., :viridis, :magma, :plasma, :cividis, :heat)

savefig("corr_hmap_6.png")