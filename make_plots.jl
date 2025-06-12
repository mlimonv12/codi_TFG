using ITensors, ITensorMPS
using Plots, Statistics
using DelimitedFiles
using CSV, HDF5
using DataFrames
gr()



function get_name(N, U, V, bdim, dt, secondorder)
    name = "converged/MPS/N_$N-U_$U-V_$V-bdim_$bdim-dt_$dt-o2_$secondorder.h5"
    return name
end



function readarray(source; unpack::Bool=false)
	df = CSV.read(source)
	if !unpack
		return convert(Matrix, df)   # returns 2-D array
	else
		return (convert(Vector, x) for x in eachcol(df))   # returns tuple of 1-D arrays
	end
end



N = 40
U = 20.0
V = 0.0
bdim = 20
dt = 0.01
secondorder = false


titol1 = get_name(N, U, V, bdim, dt, secondorder)
f = h5open(title,"r")
psi = read(f,"T",MPS)
close(f)



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



#results = readdlm("processed/CDW_vals.csv", ',', Float64)
#Vs = readdlm("processed/CDW_V.csv", ',', Float64)

#results = CSV.read("processed/CDW_vals.csv")
#Vs = CSV.read("processed/CDW_V.csv")

results = readarray("processed/CDW_vals.csv")
Vs = readarray("processed/CDW_V.csv")

println(typeof(Vs))


# CDW phase transition
Us = [4.0, 6.0, 8.0, 10.0, 15.0, 20.0]

for i in 1:6
    U = Us[i]
    if i == 1
        plot(Vs[i]/U, results[i], xlabel = "V/U", ylabel = "Δ_CDW", label = "U=$U", xlims = (0,2), lw = 3)
    else
        plot!(Vs[i]/U, results[i], xlabel = "V/U", ylabel = "Δ_CDW", label = "U=$U", xlims = (0,2), lw = 3, legendfontsize = 14)
    end
end

savefig("cdw1.png")




# CDW critical point

for i in 1:6
    nderx, ndery = num_der(Vs[i], results[i])
    
    U = round(Us[i], digits = 0)

    if i == 1
        plot(nderx, ndery/U, xlabel = "V/U", ylabel = "dΔ/dV", label = "U=$U", lw = 2)
    else
        plot!(nderx, ndery/U, xlabel = "V/U", ylabel = "dΔ/dV", label = "U=$U", lw = 2)
    end
end

savefig("cdw_numder.png")
