using ITensors, ITensorMPS

N = 4  # Unit cell size
sites = siteinds("Electron", N) # Using Qubit sites for simplicity
bdim = 5 # Bond dimension

# 2. Initialize the iMPS as a Vector{ITensor} with random tensors
links = [Index(bdim, "Link,n=$i") for i in 0:N]
A = [randomITensor(links[i], sites[i], links[i+1]) for i in 1:N]
# A = normalise_unitcell(A)



# Generates an iMPS ITensor in Vidal form with a given initial state. 1:empty, 2:up, 3:dn, 4:updn
function init_iMPS(state::Vector; bdim=16)
    
    references = Dict(0=>"Emp", 1=>"Up", 2=>"Dn", 3=>"UpDn")
    
    # Initializes iMPS with two-site unit cell 
    sites = siteinds("Electron", 2)
    links = [Index(bdim, "link-$i") for i in 0:5]
    psi = [ITensor() for _ in 1:5]

    for i in 1:5
        if iseven(i) # Γ
            T = ITensor(links[i], sites[Int(i/2)], links[i+1])
            T[links[i]=>bdim, sites[Int(i/2)]=>references[state[Int(i/2)]], links[i+1]=>bdim] = 1.0
            psi[i] = T
        else # Λ
            psi[i] = randomITensor(links[i], links[i+1])
        end
    end

    #psi[5] = psi[1] # Facillitates operations later on

    return psi, sites, links
end



# Finds the inverse of a square diagonal two-index tensor
function inv_diagonal(S::ITensor)

    Sinv = copy(S)
    s1, s2 = inds(Sinv)
    values = Int(size(S)[1])
    
    for k in 1:values
        Sinv[s1=>k, s2=>k] = 1 / (S[s1=>k, s2=>k])
    end

    return Sinv
end




# Contracts over all tensors of an iMPS to prepare for multiple observable computation. Allows to choose central Λ in case of iTEBD
function contract_full_iMPS(psi::Vector{ITensor}; central = 1)

    unitcell = copy(psi)
    
    # Find indices of all sites, for comfort
    l1_lamb, r1_lamb = inds(unitcell[1])
    l1, p1, r1 = inds(unitcell[2])
    l2_lamb, r2_lamb = inds(unitcell[3])
    l2, p2, r2 = inds(unitcell[4])

    if central == 1
        lateral_clone = unitcell[3]

        # Replace indices to allow for different iMPS contraction
        unitcell[3] = prime(unitcell[3], l2_lamb) # Avoid looping, now \Lambda_2 is the tensor on both sides
        unitcell[4] = replaceinds(unitcell[4], (r2 => l1_lamb))

        # Output tensor will have internal indices of l2_lamb', r2_lamb
        # Output lateral_clone will have indices r1, l2_lamb
    else
        lateral_clone = unitcell[1]
        
        # Prepare the lateral clone to be absorbed through the left
        lateral_clone = replaceinds(lateral_clone, (l1_lamb => r2))
        unitcell[1] = prime(unitcell[1], l1_lamb) # For consistency with previous result

        # Output tensor will have internal indices of l1_lamb', r1_lamb
        # Output lateral_clone will have indices r2, l1_lamb
    end
 

    # Contract over all tensors in the unit cell: according to this program's construction this accounts for lateral environment effects
    # (Λ_2 is duplicate in the unitcell vector, and it suffices to represent the infinite environment in Vidal form)

    # Carefully picking the order of contraction will output a tensor with the most comfortable index ordering
    # (virtual-physical-physical-virtual)
    init_contract = mod(central, 2)*2 + 1
    contracted_iMPS = unitcell[init_contract]

    for i in 1:3
        nextsite = mod1(init_contract + i, 4)
        contracted_iMPS *= unitcell[nextsite]
    end

    # Finally contract over the lateral clone at the rightmost position
    contracted_iMPS *= lateral_clone
    

    # Returns a fully internally-contracted iMPS unit cell
    return contracted_iMPS, lateral_clone
end


prova1 = randomITensor(links[1], links[2])

U, S, V = svd(prova1, links[1])

S *= 2
println(size(S))

# Per normalitzar:
#S /= norm(S)
mida = Int(size(S)[1])
S_Lind, S_Rind = inds(S)

sumatot1 = 0.0
for k in 1:mida
    sumatot1 += S[S_Lind=>k, S_Rind=>k]^2
end


S /= sqrt(sumatot1)
sumatot = 0.0

for k in 1:mida
    sumatot += S[S_Lind=>k, S_Rind=>k]^2
end
println(inds(S))