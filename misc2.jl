using ITensors, ITensorMPS
using Plots
gr()

N = 5
bdim = 10
sites = siteinds("Electron", N)
#sites = [Index(4, "link-$i") for i in 1:N]
links = [Index(bdim, "link-$i") for i in 1:N]

# Generate iMPS explicitly
psi = [randomITensor(links[mod1(i-1,N)], sites[i], links[i]) for i in 1:N]

#println(size(psi[1])[1])



# Canonicalisation via QR & normalisation test 1

j = 1

site_tensor = psi[j]
l, p, r = inds(site_tensor)

println("Left index: ", l)
println("Physical index: ", p)
println("Right index: ", r)

# Reshape D,d,D into a D*d,D

comb_lp = combiner(l, p)

site_matrix = site_tensor * comb_lp


println("Test1: ", inds(site_matrix))

# Perform QR decomposition
Q, R = qr(site_matrix, inds(site_matrix)[1])


println("Q inds: ", inds(Q))
println("R inds: ", inds(R))
println("Type of R: ", typeof(R))
prova1 = Q * R
println("Norma recuperada: ", norm(site_matrix - prova1))

#data = rand(21,100)
#Rmatrix = dense(R)
Rmatrix = Array(R, inds(R))

heatmap(1:size(Rmatrix,1),
    1:size(Rmatrix,2), Rmatrix,
    c=cgrad([:blue, :white, :red]),
    xlabel="x values", ylabel="y values",
    title="R matrix coefficients")

savefig("ViewRmatrix.png")


#prova1 = op("Cdagup", sites[j]) * op("Cup", sites[mod1(j+1, N)])