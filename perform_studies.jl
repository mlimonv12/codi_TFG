using ITensors, ITensorMPS  # Tensor elements
using Plots, Statistics     # Data processing
using HDF5, Serialization   # Data exportation
gr()

include("TEBD_EVOL.jl")

println("BEGIN STUDY")

# Modify ITensor order warning threshold
ITensors.set_warn_order(40)
#ITensors.disable_warn_order()

#sites = siteinds("Electron", N; conserve_qns=true)
#bdims = [30, 40]#[2, 4, 6, 8, 10, 12, 16, 20, 40]

# Given a fixed bond dimension, iterate over different N values to find convergence to infinite results


# An alternative approach would be calculating the variance of an observable, to see how
# Boundary conditions affect global measurements. Still, this could be mitigated by simply
# adding errorbars to measurements


#ens = [2, 6, 10, 14]#[8, 12, 16, 20, 24, 28, 32, 36]
#enedef = 
#bdims = [2, 3, 4, 5, 6, 7, 8, 10, 12, 16, 20, 25, 30, 60, 80]
#Us = [1.0, 2.0, 3.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 14.0, 16.0, 20.0, 24.0]
#Vs = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 10.0, 12.0, 14.0]
#Vs = [9.0, 11.0, 13.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0]
#[2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0]
#[1.9, 1.92, 1.94, 1.96, 1.98, 2.00, 2.02, 2.04, 2.06, 2.08, 2.1, 2.12, 2.14, 2.16]
U = 4.0

#for Vx in 1:5
    # println("Bdim: $bdim")
#    Vinst = round((Vx*0.02 + 0.6)*U, digits = 2)
#discard = [0, 5, 10, 15]

dts = [0.02]



TEBD_EVOL(
    # Number of sites
    10,
    # Hamiltonian parameters
    t = 1.0,
    U = 4.0,  
    V = 0.0,
    # Simulation parameters
    dt = 0.01,
    maxsteps = 3000,
    checkevery = 20,
    framespersecond = 10,
    cutoff = 1e-8,
    bdim = 20,
    energy_precision = 1e-30,
    secondorder = false,
    monitorobservable = false,
    savefigures = false, 
    saveconvergence = true,
    saveMPS = true
)
#    end
