println("AAA", 3.14)


U = 4.0
for Vx in 1:4
    Vinst = round((Vx*0.02 + 0.6)*U, digits = 2)
    println(Vinst)
end