using ITensors, ITensorMPS

# Computes the expected value of an operator in an iMPS Vector{ITensor} object
function imps_expect(iMPS::Vector{ITensor}, operator::ITensor, site::Int64)

    # Find number of sites the operator acts on
    s = Int64(length(size(operator))/2)
    print("\n S = ", s, "\n")

    N = length(iMPS)

    site_tensor = iMPS[site]

    # If the operator is of more than one site, absorb the (s-1) following sites into 
    # a single tensor that can be multiplied by the operator
    if s > 1
        for i in 1:(s-1)
            #print(inds(iMPS[site]), "\n")
            #print(inds(iMPS[site + mod1(i, N)]), "\n")
            site_tensor *= iMPS[mod1(site + i, N)]
        end
    end


    # Apply operator gate to sites of interest to create a tensor of dimension D,D,D,D
    #applied_tosite = iMPS[site]' * operator * iMPS[site]
    applied_tosite = site_tensor' * operator * site_tensor
    # print("cacota3: ", inds(applied_tosite), "\n\n")

    # Now a new Vector{ITensor} object is created of length N-(s-1), where s
    # is the number of sites the operator acts upon. In one of its positions
    # the object calculated above will be saved, and in the others the contractions
    # of the sites at both MPS sites (site' * site) will be stored. Finally the
    # whole iMPS will be contracted vertically
    cont_iMPS = [ITensor() for _ in 1:(N - (s - 1))]

    # We save the applied_tosite tensor to the first position of cont_iMPS
    # We don't need to preserve the position of the original iMPS, just the order.
    # cont_iMPS will store first the applied_tosite tensor and then the following
    # iMPS symmetrical sites, contracted, ordered by looping over all the iMPS from the
    # site immediately below the contracted operator up to the one immediately above the
    # site entered as an argument in this function
    cont_iMPS[1] = applied_tosite

    # For each position of this new vector that is not storing the tensor <ϕi|O|ϕi>:
    for item in 1:N-s
        # Set the cursor at position pos: the site acted upon + 1 + (s-1) (sites absorbed by multisite operator)
        pos = mod1(site + item + (s-1), N)
        pos_cont = mod1(item + 1, N - (s-1))
        indsite = "Electron,Site,n="*string(pos) # Physical index tag at given site to unprime

        # Contract symmetric iMPS sites upon which no operator acted and store at pos
        cont_iMPS[pos_cont] = noprime(iMPS[pos]', indsite) * iMPS[pos]
    end

    # Contract rest of iMPS
    expect = cont_iMPS[1] # Any position is valid: all indices are connected and will be contracted
    for i in 2:(N - (s-1))
        expect *= cont_iMPS[i]
    end

    return expect
end


# Sets a given FH state for an MPS. Input is the MPS and a vector where 1:empty, 2:up, 3:dn, 4:updn
function set_FHstate(MPS::Vector{ITensor}, sites::Vector, links::Vector, state::Vector; bdim=10)

    N = length(MPS)

    references = Dict(0=>"Emp", 1=>"Up", 2=>"Dn", 3=>"UpDn")

    for i in 1:N
        #print("AAA", typeof(state[i]), "\n")
        T = ITensor(links[i], sites[i], links[mod1(i+1,N)])
        T[links[i]=>bdim, sites[i]=>references[state[i]], links[mod1(i+1,N)]=>bdim] = 1.0
        MPS[i] = T
        #MPS[i][links(i), sites[i]=>Dict(state[i]), links(mod1(i+1,N))] = 1.0
    end

    return MPS
end


let
    N = 5
    d = 2 # Local hilbert space dimension
    χ = d^5
    sites = siteinds("Electron", N)
    #print(typeof(sites))
    #sites = [Index(4, "link-$i") for i in 1:N]
    links = [Index(10, "link-$i") for i in 1:N]

    # Generate iMPS explicitly
    #psi = [randomITensor(links[mod1(i-1,N)], sites[i], links[i]) for i in 1:N]
    #print(inds(psi[2]), "\n")

    psi = [ITensor() for _ in 1:N]

    N = 5
    sites = siteinds("Electron", N)
    #sites = [Index(4, "link-$i") for i in 1:N]
    links = [Index(10, "link-$i") for i in 1:N] # No entanglement

    for i in 1:N
        left = links[mod1(i, N)]
        right = links[mod1(i+1,N)]
        s = sites[i]
        # Initialize ITensor with bond and physical index
        T = ITensor(left, s, right)
        T[left => 10, s => "Emp", right => 10] = 1.0  # set "Emp" state explicitly
        if i==1 || i==2 || i==5# || i==5 || i==3
            T[left => 10, s => "Emp", right => 10] = 0.0  # set "Emp" state explicitly
            T[left => 10, s => "Up", right => 10] = 1.0
        end
        psi[i] = T
    end

    psi0 = copy(psi)

    psi = set_FHstate(psi, sites, links, [2,1,0,0,0])

    

    j = 1

    gate1 = op("Ntot", sites[j])
    gate2 = op("Nup", sites[j])*op("Nup", sites[mod1(j+1, N)])

    gatetot = ITensor()
    #for i in 1:N
        #gatetot ⊗= op("Ntot", sites[i])
    #end

    expect2 = scalar(imps_expect(psi, gate1, j))

    print("Escalar final = ", expect2, "\n")


    #provaa1 = op("Nup", sites[1]) * op("Nup", sites[2])
    #provaa2 = op("Nup", sites[1]) + op("Nup", sites[2])
    #provaa2 = op("Nup", sites[1]) * op("Id", sites[2:N])
    #print("\n\nprovaa1 = ", inds(provaa1), "\n")
    #print("\n\nprovaa2 = ", inds(provaa2), "\n")
    #reindex!(provaa2, [1,3,4,5,6,2,7,8,9,10])
    #provaa2 = op("Nup", sites[1]) * op("Id", sites[2:N])
    #print("\n\nprovaa1 = ", inds(provaa1), "\n")
    #print("\n\nprovaa2 = ", inds(provaa2), "\n")

    # YIPPIE!!!!!!!!!!! Contracció d'operador d'1 site feta. Cal generalitzar a 2 (i 3 hehehe)



    #prova2 = op("Cdagup", sites[1])*op("Cup", sites[mod1(1 + 1, N)])
    #print(inds(prova2), "\n")
    #prova3 = op("Nup", sites[1])
    #prova4 = op("Cdagup", sites[2])*op("Cup", sites[mod1(3, N)])# + op("Cdagup", sites[2])*op("Cup", sites[mod1(1, N)])
    #print("MIREM LA MIDA 2: ", length(size(prova2))/2, "\n")
    #print("MIREM LA MIDA 3: ", length(size(prova3))/2, "\n")
    #print("MIREM LA MIDA 4: ", length(size(prova4))/2, "\n")

    #nom1 = "operador11"

    #print("\n\n No es reconeix l'operador següent: ", nom1, ", avortant programa\n")
    #print(length(sites))
    #print(typeof(prova2))
    #throw(InterruptException())
    #scalar = 2
    #print(scalar)

end