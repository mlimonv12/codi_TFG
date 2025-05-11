using ITensors, ITensorMPS

println("Funciona")

N = 5

sites = siteinds("Electron", N)

j = 1

mirar1 = op("Cdagdn", sites[j])

t = 1.0

os = OpSum()
os += -t, "Cdagup", 1, "Cup", 2
os += -t, "Cup", 1, "Cdagup", 2

H = MPO(os, sites)

println(mirar1)


#println("\n", H[1])