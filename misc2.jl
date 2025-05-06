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

j = 2

prova1 = op("Cdagup", sites[j]) * op("Cup", sites[mod1(j+1, N)])

println(inds(prova1))