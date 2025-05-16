# Generates an iMPS ITensor in Vidal form with a given initial state. 1:empty, 2:up, 3:dn, 4:updn
function init_iMPS(state::Vector; bdim=16)
    
    references = Dict(0=>"Emp", 1=>"Up", 2=>"Dn", 3=>"UpDn")
    
    # Initializes iMPS with two-site unit cell 
    sites = siteinds("Electron", 2)
    links = [Index(bdim, "link-$i") for i in 1:4]
    psi = [ITensor() for _ in 1:4]

    for i in 1:4
        if iseven(i) # Γ
            T = ITensor(links[i], sites[Int(i/2)], links[mod1(i+1, 4)])
            T[links[i]=>1, sites[Int(i/2)]=>references[state[Int(i/2)]], links[mod1(i+1, 4)]=>1] = 1.0
            psi[i] = T
            #psi[i] = randomITensor(links[i], sites[Int(i/2)], links[i+1])

        else # Λ
            psi[i] = ITensor(zeros(bdim, bdim), links[i], links[mod1(i+1, 4)])
            for k in 1:bdim
                psi[i][links[i]=>k, links[mod1(i+1, 4)]=>k] = 1.0
            end
        end
    end

    psi = normalise_unitcell(psi)

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


# Normalises a square diagonal two-index tensor
function normalise_diagonal(S::ITensor)

    mida = Int(size(S)[1])
    S_Lind, S_Rind = inds(S)

    sumatot = 0.0
    for k in 1:mida
        sumatot += S[S_Lind=>k, S_Rind=>k]^2
    end

    return S / sqrt(sumatot)
end

# Normalises Λ tensors of an iMPS unit cell
function normalise_unitcell(unitcell::Vector{ITensor})

    bdim = Int(size(unitcell[1])[1])

    unitcell[1] = normalise_diagonal(unitcell[1])*sqrt(bdim)
    unitcell[3] = normalise_diagonal(unitcell[3])*sqrt(bdim)
    #unitcell[2] /= norm(unitcell[2])
    #unitcell[4] /= norm(unitcell[4])
    # unitcell[5] = unitcell[1]

    return unitcell
end


# Finds the site index of a Γ tensor
function find_siteindex(sitetensor::ITensor)

    for index in inds(sitetensor)
        if hastags(index, "Site")
            return index
        end
    end
    return
end


# Contracts over all tensors of an iMPS to prepare for multiple observable computation. Allows to choose central Λ in case of iTEBD
function contract_full_iMPS(psi::Vector{ITensor}; central = 1, matrix = false)

    unitcell = copy(psi)
    # New attempt: the link indices are defined in a cyclic way, spares contraction conditions

    # Central == 1 means that the matrix Λ_1 is taken to be the centre of orthogonality of this contraction

    if central == 1     # Λ_2 - Γ_2 - Λ_1 - Γ_1 - Λ_2
        rightmost_clone = copy(unitcell[3])

        unitcell[3] = prime(unitcell[3], commonind(unitcell[2], unitcell[3])) # Leftmost contraction
        rightmost_clone = prime(rightmost_clone, commonind(rightmost_clone, unitcell[4]))

        # Output tensor will have internal indices of l2_lamb', r2_lamb'

        # Find combiners for transfer matrix form
        if matrix
            left_internal = prime(commonind(unitcell[2], rightmost_clone)) # Primed, as we have just done
            left_physical = find_siteindex(unitcell[4])
            right_internal = prime(commonind(unitcell[3], unitcell[4]))
            right_physical = find_siteindex(unitcell[2])

            cleft = combiner(left_internal, left_physical)
            cright = combiner(right_internal, right_physical)
        end

    else                # Λ_1 - Γ_1 - Λ_2 - Γ_2 - Λ_1
        rightmost_clone = unitcell[1]
        
        unitcell[1] = prime(unitcell[1], commonind(unitcell[1], unitcell[4]))
        rightmost_clone = prime(rightmost_clone, commonind(rightmost_clone, unitcell[2]))

        # Output tensor will have internal indices of l1_lamb', r1_lamb'

        # Find combiners for transfer matrix form
        if matrix
            left_internal = prime(commonind(unitcell[4], rightmost_clone)) # Primed, as we have just done
            left_physical = find_siteindex(unitcell[2])
            right_internal = prime(commonind(unitcell[1], unitcell[2]))
            right_physical = find_siteindex(unitcell[4])

            cleft = combiner(left_internal, left_physical)
            cright = combiner(right_internal, right_physical)
        end
    end
 

    # Contract over all tensors in the unit cell: according to this program's construction this accounts for lateral environment effects
    # (Λ_2 is duplicate in the unitcell vector, and it suffices to represent the infinite environment in Vidal form)
    contracted_iMPS = unitcell[1]

    for i in 2:4
        contracted_iMPS *= unitcell[i]
    end

    # Finally contract over the lateral clone at the rightmost position
    contracted_iMPS *= rightmost_clone
    
    # Prime both indices of rightmost_clone to facillitate inversion later on
    noprime!(rightmost_clone)
    prime!(rightmost_clone)

    if matrix
        return contracted_iMPS * cleft * cright, rightmost_clone, cleft, cright
    end

    # Returns a fully internally-contracted iMPS unit cell
    return contracted_iMPS, rightmost_clone
end

