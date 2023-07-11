using SparseArrays, LinearAlgebra, Random, MAT, PyPlot, CSV, Tables, CavityTools
#Random.seed!(1234)


function tempo()
    # parameters
    densit = []
    N = 10000                  #number of nodes (a square of a number for the square lattice)                      
    ncol = 1                      #number of protein types (colors)
    type_lattice = "square"         #type of lattice (ring, square, triangle, hexagon, moore)
    T = 1e7                         #final time
    g = 50                       #aggregation constant (g >1 homologous neighbours reduce the prob )
    kI = 1e-6                     #insertion coefficient
    kE = 1e20                      #extraction coefficient 
    thrE = 50 #trunc(Int64, (pi/2)*9) #extraction threshold
    soglia = 20

    t = 0.                        #initial time
    rho = 0.0                    #constant density

    function create_lattice(type_lattice, N, kD=1)
     """create_lattice generates different types of lattice: ring, square, triangle"""  


        if type_lattice == "ring"
            V = 2
            nodes = collect(1:N)
            nodes =reshape(nodes, N, 1)
            left_nn = circshift(nodes, 1)
            right_nn = circshift(nodes, N-1)
            nn = hcat(left_nn,right_nn)


            v = ones(Int64, N)*V
            kD = 1 #never used
            return nn, v, V, N, kD
        end


        if type_lattice == "square"
            V = 4
            L = sqrt(N)
            nn = zeros(Int64, N, V)
            @inbounds for i in 1:N
                nn[i,1] = mod(i-2,L)+1 +((i-1)÷L)*L    #left 
                nn[i,2] = mod(i,L) +1 + ((i-1)÷L)*L    #right
                nn[i,3] = mod(i-1-L,N) +1              #up
                nn[i,4] = mod(i-1+L,N) +1              #down
            end 

            v = ones(Int64, N)*V
            kD = 1 
            return nn, v, V, N, kD
        end


        if type_lattice == "triangle"
            l4 = trunc(Int64, sqrt(N))   #square with side l4

            c3=2*round(Int64, 2^(-1.)*3^(1/4.)*l4)
            r3=round(Int64, 3^(-1/4.)*l4)
            N3=c3*r3 
            println("col: " ,c3, "   row: ", r3, "    nodes: ", N3)

            nodes = collect(Int64, 1:N3)
            ind = reshape(nodes, r3, c3)

            r = collect(1:r3)       
            c = collect(1:c3)   

            codd = c[(c.%2).==1]
            ceven = c[(c.%2).==0]

            V = 3  #num of neighbours cells
            nn = zeros(Int64, N3, V)

            ### #down triangle (odd col)
            for o in ind[:, codd] 
                nn[o,1] = (o + r3 - 1) + r3*(mod(o-r3-2,r3)÷(r3-1))          #up
                nn[o,2] = o + r3                                             #right-down 
                nn[o,3] = mod(o-1 - r3,r3*c3) +1                             #left-down                                                   #U
                #println("o: ", o ,"      nn[0,1]: ", nn[o,1],"      nn[0,2]: ", nn[o,2],"      nn[0,3]: ", nn[o,3]) 
            end    

            #up triangle (even col)
            for o in ind[:, ceven]   
                nn[o,1] = o - r3                                            #left-up 
                nn[o,2] = mod(o + r3,r3*c3)                                 #right-up
                nn[o,3] = (o - r3 + 1) - r3*(mod(o-r3-1,r3)÷(r3-1))         #down                
                #println("o: ", o ,"      nn[0,1]: ", nn[o,1],"      nn[0,2]: ", nn[o,2],"      nn[0,3]: ", nn[o,3])
            end    

            v = ones(Int64, N3)*V
            kD = 3*sqrt(3)*kD/4 
            return nn, v, V, N3, kD
        end


        if type_lattice == "hexagon"
            l4=trunc(Int64, sqrt(N))   #square with side l4
            c6=2*round(Int64, 2^(-1/2.)*3^(-1/4.)*l4)
            r6=round(Int64, 2^(-1/2.)*3^(1/4.)*l4)
            N6=c6*r6 
            println("col: " ,c6, "   row: ", r6, "    nodes: ", N6)

            nodes = collect(Int64, 1:N6)
            ind = reshape(nodes, r6, c6)

            r = collect(1:r6)       
            c = collect(1:c6)   

            codd = c[(c.%2).==1]
            ceven = c[(c.%2).==0]

            V = 6  #num of neighbours cells
            nn = zeros(Int64, N6, V)

            ### #down triangle (odd col)
            for o in ind[:, codd] 
                nn[o,1] = mod(o - r6 -1,r6*c6 ) +1                            #ul
                nn[o,2] = mod(o - 1, r6*c6 ) +r6*(mod(o-r6-2,r6)÷(r6-1))      #u 
                nn[o,3] = o + r6                                              #ur   
                nn[o,4] = o + r6 +1  - r6*(mod(o-r6-1,r6)÷(r6-1))             #dr
                nn[o,5] = mod(o + 1, r6*c6 ) - r6*(mod(o-r6-1,r6)÷(r6-1))     #d
                nn[o,6] = mod(o -1 -r6,r6*c6 ) +2 - r6*(mod(o-r6-1,r6)÷(r6-1))     #dl
                #println("o: ", o ,"      nn[0,1]: ", nn[o,1],"      nn[0,2]: ", nn[o,2],"      nn[0,3]: ", nn[o,3]) 
            end    

            #up triangle (even col)
            for o in ind[:, ceven]   
                nn[o,1] = mod(o - r6 -1,r6*c6+1 ) + r6*(mod(o-r6-2,r6)÷(r6-1))    #ul 
                nn[o,2] = mod(o - 1, r6*c6 ) +r6*(mod(o-r6-2,r6)÷(r6-1))          #u 
                nn[o,3] = mod(o +r6  -1 , r6*c6 ) +r6*(mod(o-r6-2,r6)÷(r6-1))     #ur
                nn[o,4] = mod(o + r6, r6*c6)                                      #dr
                nn[o,5] = mod(o + 1, r6*c6 +2) - r6*(mod(o-r6-1,r6)÷(r6-1))       #d
                nn[o,6] = o - r6                                                  #dl
                #println("o: ", o ,"      nn[0,1]: ", nn[o,1],"      nn[0,2]: ", nn[o,2],"      nn[0,3]: ", nn[o,3]) 
            end 

            v = ones(Int64, N6)*V
            kD = sqrt(3)*kD/2 
            return nn, v, V, N6, kD 
        end 

        if type_lattice == "moore"
            V = 8
            L = sqrt(N)
            nn = zeros(Int64, N, V)
            nodes = collect(1:N) 
            perm_L = circshift(nodes, L)
            perm_Lm1L = circshift(nodes, (L-1)*L)
            for i in 1:N
                nn[i,1] = mod(perm_L[i]-2,L)+ 1 +((perm_L[i]-1)÷L)*L #up-left
                nn[i,2] = mod(i-1-L,N) +1                            #up
                nn[i,3] = mod(perm_L[i],L) +1 + ((perm_L[i]-1)÷L)*L  #up-right   
                nn[i,4] = mod(i,L) +1 + ((i-1)÷L)*L                  #right
                nn[i,5] = mod(perm_Lm1L[i],L) +1 + ((perm_Lm1L[i]-1)÷L)*L    #down-right
                nn[i,6] = mod(i-1+L,N) +1                                    #down
                nn[i,7] = mod(perm_Lm1L[i]-2,L)+1 +((perm_Lm1L[i]-1)÷L)*L    #down-left 
                nn[i,8] = mod(i-2,L)+1 +((i-1)÷L)*L                          #left 
            end

            v = ones(Int64, N)*V
            kD = 2*kD/3 
            return nn, v, V, N, kD  
        end

    end

    function rate(mp,g::Int64=g, kD=kD)           #mp = moves possible
        return kD*mp[2]/g^mp[1]
    end 
    
    nn, v, V, N, kD = create_lattice(type_lattice, N);
    mp = [[h,e] for h=0:V+1 for e=0:V+1 if ((h+e < V+1) & (e>0))];
    mpi = Dict(i => c for (c,i)=enumerate(mp))
    rates = map(rate, mp);
    
    eta = rand(1:ncol, N).*(rand(N).<rho);
    
    #density 
    function counting(eta)
        nums = Int64[]
        for i in sort!(unique(eta))
            n = count(x->x==i,eta)
            push!(nums, n)
        end
        return nums
    end

    nums = counting(eta) #count the number of particles for each species in starting from 0
    if length(nums) == 1
        rho_tot = 0.
    else
        rho_tot = sum(nums[2:end])/N
    end
    #println("nums " ,nums)
    rho_i = zeros(Float64, 1, ncol)

    for (ind, col) in enumerate(nums[2:end])
        #println(ind,"  ", col)
        rho_i[ind] = nums[2:end][ind]/N
    end
    rho_tot, rho_i
    
    ne = (v  - sum(eta[nn].>0, dims=2)).*(eta.>0)
    eta1 = repeat(reshape(eta, N ,1), 1, V)
    nh = sum((eta1 .== eta[nn]) .& (eta1.>0) ,dims = 2);  #array of homologous neighbours for each node
    #println("number of empty neighbors ", ne, "\n")
    #println("eta1 is a multi-column-version of the eta matrix: \n", eta1, "\n")
    #println("number of homologhe neighbors ", nh)

    m = hcat(nh,ne);
    
    sd = Vector{Int64}[]
    isd = -ones(Int64, N)
    @inbounds for i in 1:length(mp)
        push!(sd, Int64[])
        for j in 1:N           # N = length(m)
            if all(m[j,:] .== mp[i])
                append!(sd[i],j) 
                isd[j] = length(sd[i])
            end
        end
    end

    #nodes where threre is not particle prensent
    pel = -ones(Int64, N)                       #pel = position of the empty nodes in the empty list
    el = Int64[]
    for e in 1:N           #This loop can be parallelized                
        if eta[e] == 0
            push!(el, e)
            pel[e] = length(el)
        end
    end

    nₚ = 0.   # number of monomers
    Nd = 0    # number of domains
    nd = 0.   # number of molecules inside the clusters (density)
    ns = 0.   # number of molecules inside the surface of a clusters (density)
    Nin = 0   #number of insertions
    Next = 0  # number of extracted domains
    Ndis = 0  # number of evaporations (domain that evaporates)
    next = 0. # number  of molecules inside an extracted cluster (is approximately equal to Next*thrE)
    τₚ = 0.    #productive domain timelife
    τᵤ = 0.    #unproductive domain timelife
    pos = 0
    lcc = Vector{Int64}[]        #= list of connected components without distinction of the col.
                    Each sub-array contains the position of the particles belonging to that cluster into i =#
    cc = -ones(Int64, N)     # it indicates the cluster (sub-array) where it is the particle
    icc = -ones(Int64, N)    # it gives us the position of the particle inside the sub-arrays 

    lccS = Vector{Int64}[] #perimeters of the connected components
    iccS = -ones(Int64, N) # it gives us the position of the particle inside the sub-arrays
    tin = Float64[]
    for i in 1:N
        if (eta[i] != 0) & (cc[i] == -1)
            #println(i)
            pos += 1
            cc[i] = pos
            icc[i] = 1
            push!(lcc, [i])
            l = length(lcc)
            if sum(eta[nn[i,:]] .== eta[i]) != V  #particle on the surface
                iccS[i] = 1
                push!(lccS, [i])
            else
                push!(lccS, [])
            end
            @inbounds for j in lcc[l]
                @inbounds for k in nn[j,:]
                    if (eta[k] == eta[i]) & (cc[k] == -1)
                        #print(nn[j,:])
                        cc[k] = pos
                        push!(lcc[l], k)
                        icc[k] = length(lcc[l]) 

                        if sum(eta[nn[k,:]] .== eta[i]) != V  #particle on the surface
                        push!(lccS[l], k)
                        iccS[k] = length(lccS[l])
                        end

                    end  
                end
            end
            if length(lcc[l]) > 1
                push!(tin, t)
                Nd += 1
                nd += length(lcc[l])/N  #num of particles inside a cluster
                ns += length(lccS[l])/N #number of particles in the surface of a cluster
            else
                nₚ += 1/N
                push!(tin, -1.)
            end
        end    
    end

    #clusters over size
    over_thrE = findall(x->length(x)>=thrE, lcc);
    
    function choose_move(qs)
        rd = rand(Float64)
        move = 1
        #println(rd)
        while move <= length(qs)
            #println(move)
            if rd < qs[move]
                if move < length(qs) - 1
                    return move, "diffusion"
                elseif move < length(qs)
                    return move, "insertion"
                else 
                    return move, "extraction"
                end    
            else
                move += 1
            end    
        end
    end
    
    function diffusion(move, sd, eta, nn=nn::Array{Int64,2}) 
        pc = rand(sd[move])                                  # 2) random choice of the particle to move
        rc = rand(nn[pc,:][findall(x->x==0, eta[nn[pc,:]])]) # 3) select an empty neighbours of pc at random
        #update the system
        eta[rc] = eta[pc]
        eta[pc] = 0

        return pc, rc, eta
    end
    
    function insertion(el, eta, ncol::Int64=ncol)
        nc = rand(el)          # 2) choose a node at random (nc = chosen node)
        col = rand(1:ncol)      # 3) choose the protein type to insert
        #update the system
        eta[nc] = col

        return nc, col, eta
    end
    
    function extraction(over_thrE, eta, lcc=lcc::Array{Array{Int64,1},1})
        ceot = rand(over_thrE)         #ceot = select a cluster from over_thrE (position of the cluster in sd)

        col = eta[lcc[ceot]][1]        #color of the cluster that has to be extracted 

        #update the system
        eta[lcc[ceot]] .= 0  #extraction

        return ceot, col, eta
    end 
    
    function create_m_loc(lnm, eta, nn=nn::Array{Int64,2}, v=v::Array{Int64,1}, V=V::Int64)
        #common part
        ne_loc = (v[lnm]  - sum(eta[nn][lnm,:].>0, dims=2)).*(eta[lnm].>0) #array of empty neighbours for each node
        eta1 = repeat(reshape(eta[lnm], length(lnm) ,1), 1, V)
        nh_loc = sum((eta1 .== eta[nn[lnm,:]]) .& (eta1.>0) ,dims = 2);  #array of homologous neighbours for each node
        m_loc = hcat(nh_loc,ne_loc)
        #println("m_loc: \n", m_loc)
        return m_loc
    end
    
    function update_local_sd_isd(lnm, m, m_loc, sd, isd, mp=mp, mpi=mpi )
        @inbounds for i in lnm
            if m[i,:] in mp                 #this use m and not m_loc, for retrieve the particles
                ind_sd = mpi[m[i,:]]
                pos_in_sd = isd[i]   #position in sd
                sd[ind_sd][pos_in_sd] = sd[ind_sd][end]
                isd[sd[ind_sd][end]] = pos_in_sd
                isd[i] = -1
                pop!(sd[ind_sd])
            end
        end
        @inbounds for (i,x) in enumerate(lnm)
            if m_loc[i,:] in mp
                ind_sd = mpi[m_loc[i,:]]
                append!(sd[ind_sd], x) 
                isd[x] = length(sd[ind_sd])  #c'era i
            end
        end
        return sd, isd
    end 
    
    function update_local_el_pel_insertion(nc, el, pel)
        el[pel[nc]] = el[end]
        l = pop!(el)
        pel[l] = pel[nc]
        pel[nc] = -1
        return el, pel
    end
    
    function update_local_el_pel_diffusion(pc, rc, el, pel)
        el[findall(x->x==rc, el)] .= pc
        pel[pc] = pel[rc]
        pel[rc] = -1
        return el, pel
    end
    
    #=function update_local_el_pel_extraction(ceot, el, pel, lcc)
        @inbounds for i in lcc[ceot]
            push!(el, i)
            pel[i] = length(el)
        end
        return el, pel
end =#
    
    function update_local_lcc_cc_icc_insertion(lcc, cc, icc, over_thrE, over_soglia, dimeri, lccS, iccS, eta, nc, nₚ, nd, Nd, ns, tin, t, nn=nn::Array{Int64,2}, thrE::Int64=thrE, N::Int64=N)
        if eta[nc] ∉ eta[nn[nc,:]] #1) isolated particle (monomer) 
            #println("homologous neighbors = 0")
            push!(lcc, [nc])
            push!(tin, -1)
            cc[nc] = length(lcc)
            icc[nc] = 1
            push!(lccS, [nc])
            iccS[nc] = 1
            nₚ += 1/N 

        elseif sum(eta[nc] .== eta[nn[nc,:]]) == 1     # 2) homologous particle in the neighbourhood = 1
            #println("homologous neighbors = 1")
            ## I want to find the neighbour cluster
            plcc = cc[nn[nc,:]][eta[nn[nc,:]] .== eta[nc]][1]  #position in lcc
            ccn = copy(lcc[plcc])  #ccn = connected component neighbours
            #print("ccn: ", ccn)
            append!(lcc[plcc], nc)
            icc[nc] = length(lcc[plcc])
            cc[nc] = plcc
            if length(ccn) == 1 #join of two monomers
                #println("cluster: ", lcc[plcc])
                tin[plcc] = t    #the time where a new domain is formed
            end       

            an = nn[nc,eta[nn[nc,:]] .== eta[nc]][1]                              #  
            #println("adjecient node to nc (an): ", an)                       #  
            if sum(eta[nn[an,:]] .== eta[nc]) == V  #particle in the interior  #
                lccS[cc[an]][iccS[an]] = nc                                   #   
                iccS[nc] = iccS[an]                                           #
                iccS[an] = -1                                                 #
            else                                                              #    
                append!(lccS[plcc], nc)                                       #  perimeter 
                iccS[nc] = length(lccS[plcc])                                 # 
            end                                                               #  

            #update over_thrE
            if (length(lcc[plcc]) >= thrE) && (plcc ∉ over_thrE)
                push!(over_thrE, plcc)
            end
            
            if length(ccn) == 1
                nₚ -= 1/N
                Nd += 1
                nd += length(lcc[plcc])/N
                ns += length(lccS[plcc])/N  #forse è più efficiente con 2/N
            else
                nd += 1/N
                if sum(eta[nn[an,:]]) != V  
                    ns += 1/N                
                end
            end
            #print("np: ", nₚ, "   nd: ", nd, "  Nd: ", Nd)

        else
            #3) the new particle can marge two or more clusters (nh[nc] >= 2)
            #println("homologous neighbors => 2")
            cn = cc[nn[nc,:][eta[nn[nc,:]] .== eta[nc]]] #which clusters belonging to the homologous neighbors
            #println("lcc[cn]: ", lcc[cn])
            clusters = unique(cn)

            len_cl = map(x -> length(x), lcc[clusters]) 
            num_dom = sum(len_cl .> 1)
            #print("num_dom", num_dom)

             #different clusters in the neighborhood 
            clusters = sort!(clusters, rev=true)
            #println("clusters: ", clusters)

            #println(clusters)
            merged_cl = [nc]  #the new bigger cluster has to contain the particle nc
            merged_clS = []
            if sum(eta[nn[nc,:]] .== eta[nc]) != V
                push!(merged_clS, nc)
            end
            times_rm = Float64[]
            @inbounds for cl in clusters
                if length(lcc[cl]) >= thrE
                    over_thrE[cl] = over_thrE[end]
                    pop!(over_thrE)
                end
                if length(lcc[cl]) > 1      #update nd, nₚ 
                    nd -= length(lcc[cl])/N
                    ns -= length(lccS[cl])/N
                else
                    nₚ -= 1/N
                end    

                merged_cl = append!(merged_cl, lcc[cl]) # create a bigger cluster
                merged_clS = append!(merged_clS, lccS[cl])
                #print("merged_cl: ",merged_cl)

                if lcc[cl] != lcc[end]
                    cc[lcc[end]] .= cc[lcc[cl]][1]
                    cc[lcc[cl]][1] = -1
                    lcc[cl] = lcc[end]
                    tin[cl] = tin[end]
                    pop!(lcc)
                    t_rm = pop!(tin)                  #remove formation time


                    lccS[cl] = lccS[end]
                    pop!(lccS)
                else
                    cc[lcc[end]][1] = -1
                    pop!(lcc)                #remove the old cluster 
                    pop!(lccS)
                    t_rm = pop!(tin)                #remove formation time
                end
                if t_rm != -1. 
                    append!(times_rm, t_rm)
                end

            end
            #println("removed times ", times_rm)
            if times_rm == []
                t_form = t
            else
                t_form = minimum(times_rm)
            end

            #find particle in the cluster interior
            #merged_clS = copy(merged_cl)
            p_in = nn[nc,:] .* (sum(eta[nn[nn[nc,:],:]], dims=2) .== V)
            p_in = p_in[p_in .> 0]
            #println("internal particels: ", p_in)

            iccS[p_in] .= -1
            for i in p_in
                filter!(x->x≠i, merged_clS)
            end
            #println("lccS: ", lccS)

            push!(lccS, merged_clS)
            push!(lcc, merged_cl)    #add the new cluster
            push!(tin, t_form)            # formation time
            nd += length(merged_cl)/N
            ns += length(merged_clS)/N
            if num_dom == 0 
                Nd += 1
            elseif num_dom > 1
                Nd -= (num_dom -1) 
            end
            if (length(merged_cl) >= thrE) 
                push!(over_thrE, length(lcc))
            end
            
            cc[lcc[end]] .= length(lcc)
            #println(merged_cl)

            for (ind, p) in enumerate(merged_cl)
                icc[p] = ind
            end
            for (ind, p) in enumerate(merged_clS)   #
                iccS[p] = ind                       # perimeter
            end                                     #
        end
        return lcc, cc ,icc, over_thrE, over_soglia, dimeri, lccS, iccS, nₚ, nd, Nd, ns, tin
    end
    
    function update_local_lcc_cc_icc_diffusion(lnm, lcc, cc, icc, over_thrE, over_soglia, dimeri, lccS, iccS, eta, pc, rc, nₚ, nd, Nd, ns, tin, t, τᵤ, Ndis, nn=nn::Array{Int64,2}, thrE::Int64=thrE)
        #remove changed clusters
        lpcr =Int64[]           #list with the positions of the clusters to remove
        @inbounds for i in lnm
            if ((cc[i] != -1) & (eta[i] == eta[rc])) | (i == pc)  #cc[i] != -1 is equivalent to (eta[i] != 0) 
                append!(lpcr, cc[i])
            end
        end 
        lpcr = sort!(unique(lpcr), rev=true)
        #println("lpcr: ", lpcr)
        times_rm = Float64[]
        lcc_rm = []
        cl_rm = []
        @inbounds for i in lpcr
            if length(lcc[i]) >= thrE
                p = findall(x->x==i, over_thrE)[1] #update over_thrE
                over_thrE[p] = over_thrE[end]
                pop!(over_thrE)
            end
            if length(lcc[i]) > 1    #update nₚ, nd, Nd
                nd -= length(lcc[i])/N
                ns -= length(lccS[i])/N
                Nd -= 1
            else
                nₚ -= length(lcc[i])/N
            end
            if lcc[i] != lcc[end]
                cl_rm = lcc[i]  #cluster that will be removed
                t_rm = tin[i]
                tp_cc = cc[lcc[i]][1]
                cc[lcc[i]] .= -1
                cc[lcc[end]] .= tp_cc
                lcc[i] = lcc[end]
                lccS[i] = lccS[end]
                tin[i] = tin[end]
                pop!(lcc)
                pop!(lccS)
                pop!(tin)   #remove formation time
            else
                cc[lcc[end]] .= -1
                cl_rm = pop!(lcc)  #rm_cl = cluster removed
                pop!(lccS)
                t_rm = pop!(tin)
            end 
            if (t_rm != -1.) #&& (pc ∉ cl_rm)
                append!(times_rm, t_rm)
                append!(lcc_rm, [cl_rm])
            end
        end
        #println("removed times ", times_rm)
        #println("removed clusters ", lcc_rm)

        icc[pc] = -1
        iccS[nn[rc,:]] .= -1 
        #println("nₚ: ", nₚ)
        #println("nd: ", nd)
        #println("Nd: ", Nd)


        #add new clasters
        tlᵤ = Float64[]
        colors = Int64[]
        starting_points = [k for k in nn[pc, :] if eta[rc] == eta[k]]
        @inbounds for s in starting_points
            if cc[s] ∉ colors
                c = length(lcc) + 1
                cc[s] = c
                icc[s] = 1
                push!(colors, c)
                push!(lcc, [s])          
                if sum(eta[nn[s,:]] .== eta[s]) != V  #particle on the surface            
                    push!(lccS, [s])
                    iccS[s] = 1
                end
                @inbounds for j in lcc[c]
                    @inbounds for k in nn[j,:]
                        if (eta[k] == eta[j]) & (cc[k] ∉ colors)
                            push!(lcc[c], k)
                            if sum(eta[nn[k,:]] .== eta[k]) != V
                                push!(lccS[c], k)
                                iccS[k] =length(lccS[c])
                            end
                            cc[k] = c
                            icc[k] =length(lcc[c]) 
                        end
                    end               
                end
                if (length(lcc[c]) >= thrE) && (length(lcc) ∉ over_thrE) && all(x->x<=5, [length(lcc[c] ∩ over_soglia[i]) for i = 1:length(over_soglia)])        
                    push!(over_thrE, length(lcc))
                    push!(over_soglia, lcc[c])
                end 
                
                if length(lcc[c]) > 1
                    if sum([(cl ∩ lcc[c]) ⊆ lcc[c] for cl in lcc_rm if (length(cl) > 1 && !isempty(cl ∩ lcc[c]))] .> 0) > 0

                       pos = findall([(cl ∩ lcc[c]) ⊆ lcc[c] && (length(cl) > 1 && !isempty(cl ∩ lcc[c])) for cl in lcc_rm ])
                        #println("pos ", pos)
                        #println("cluster: ", lcc[c], "     time: ", minimum(times_rm[pos]), "     rm clusters: ", lcc_rm[pos])
                        push!(tin, minimum(times_rm[pos])) #add domain formation                 
                    else
                        push!(tin, t)
                        #println(lcc[c]) 
                        #println("time: ", t)
                    end
                    Nd += 1
                    nd += length(lcc[c])/N  #num of particles inside domains
                    ns += length(lccS[c])/N            
                    #println("Nd: ", Nd)
                else
                    nₚ += 1/N
                    push!(tin, -1.)             
                    # unproductive domains timelife 
                    #println("lcc[c]: ", lcc[c])
                    if length(times_rm[findall(lcc[c] .∈ lcc_rm)]) > 0 

                        τᵤ = t - times_rm[findall(lcc[c] .∈ lcc_rm)][1]
                        #println("timelife of unproductive domain: ", τᵤ)
                        push!(tlᵤ, τᵤ)

                    end     
                end
            end
        end 
        if length(tlᵤ) > 0
            τᵤ = unique!(tlᵤ)[1]                                        
            Ndis += 1
        else
            τᵤ = 0.
        end                                          

        return lcc, cc, icc, over_thrE, over_soglia, dimeri, lccS, iccS, nₚ, nd, Nd, ns, tin, τᵤ, Ndis            
    end
    
    #= function update_local_lcc_cc_icc_extraction(ceot, lcc, cc, icc, over_thrE, lccS, iccS, nₚ, nd, Nd, ns, Next, next, tin, t, τₚ)

        nd -= length(lcc[ceot])/N #update nd
        ns -= length(lccS[ceot])/N
        Nd -= 1
        Next += 1 
        next += length(lcc[ceot])/N

        if lcc[ceot] != lcc[end]
            cc[lcc[ceot]] .= -1
            icc[lcc[ceot]] .= -1
            lcc[ceot] = lcc[end]
            t_rm = tin[ceot]
            tin[ceot] = tin[end]
            cc[lcc[end]] .= ceot
            pop!(lcc)
            pop!(tin)

            iccS[lccS[ceot]] .= -1        #
            lccS[ceot] = lccS[end]        # perimeter
            pop!(lccS)                    #  
        else
            last = pop!(lcc)
            t_rm = pop!(tin)                  #remove formation time
            cc[last] .= -1
            icc[last] .= -1 

            lastS = pop!(lccS)            # perimeter
            iccS[lastS] .= -1              #
        end  
        τₚ= t - t_rm    #productive residence time

        pceot = findall(x->x==ceot, over_thrE)#position of ceot inside over_thrE
        #println("pceot: ",pceot)
        over_thrE[pceot] .= over_thrE[end] 
        pop!(over_thrE)


        #println("ok2")

        return lcc, cc, icc, over_thrE, lccS, iccS, nₚ, nd, Nd, ns, Next, next, tin, τₚ
end =#
    
    function statistics_1species(eta,t, rho_tot, nₚ, nd, Nd, ns, Next, next, cont, dt, sum_dt, sum_rho_dt, 
                                 sum_rho2_dt, Ndis, Nin, T=T, N=N, thrE=thrE, g=g, ncol=ncol, kI=kI; points = 1000000)

        sum_dt += dt
        sum_rho_dt += dt*rho_tot
        sum_rho2_dt += dt*rho_tot^2
        cont += 1

        if cont == points
            avg_rho = sum_rho_dt/sum_dt
            avg_rho2 = sum_rho2_dt/sum_dt
            var = avg_rho2 - avg_rho^2
            
            write_on_file(eta, t, avg_rho, var, Next, Nin, Ndis, n, densit, over_soglia, dimeri)
            cont = 0
            sum_dt = 0.
            sum_rho_dt = 0. 
            sum_rho2_dt = 0.
            #Next = 0
            #Ndis = 0
            #Nin = 0
        end
        return cont, sum_dt, sum_rho_dt, sum_rho2_dt, Next, Ndis, Nin
    end
    
    function write_on_file(eta, t, avg_rho, var, Next, Nin, Ndis, n, densit, over_soglia, dimeri, T::Float64=T, N::Int64=N, thrE::Int64=thrE, g::Int64=g, ncol::Int64=ncol, kI::Float64=kI)
        #data series
        
        mol_sorting = open("Mol_sorting_2D_noextr(nodes_$(N)_g_$(g)_kI_$(kI)_thrE_$(thrE))", "a")
        data = "$t" *"    " * "$avg_rho" * "    " * "$var" * "    " * "$Next" * "    " * "$Nin" * "\n"
        write(mol_sorting, data)
        close(mol_sorting)
        
        corr = open("Correlation_2D_noextr(nodes_$(N)_g_$(g)_kI_$(kI)_thrE_$(thrE)_soglia_$(soglia))", "w")
        data1 = "$over_soglia" * "\n"
        write(corr, data1)
        close(corr)
        
        #dime = open("Dimeri_2D(nodes_$(N)_g_$(g)_kI_$(kI)_thrE_$(thrE)", "w")
        #data2 = "$dim" * "\n"
        #write(dime, data2)
        #close(dime)
        
        #extr = open("Extraction_2D20(nodes_$(N)_g_$(g)_kI_$(kI)_thrE_$(thrE))", "w")
        #write(extr, "$n \n")
        #close(extr)
        
        #density = open("Density(nodes_$(N)_g_$(g)_kI_$(kI)_thrE_$(thrE))", "a")
        #write(density, "$avg_rho \n")
        #close(density)
        
        #imshow
        #system = open("imshow(ncol_$(ncol)_g_$(g)_kI_$(kI)_thrE_$(thrE)type_lattice_$(type_lattice))", "a")
        #system = "imshow(ncol_$(ncol)_g_$(g)_kI_$(kI)_thrE_$(thrE)type_lattice_$(type_lattice)).jld"
        #save(system, eta, "w")
        #write(system, "$eta \n" )
        #close(system)
    end
    
    #statistics
    cont = 0
    sum_dt = 0.          
    sum_rho_dt = 0.       #rho*t
    sum_rho2_dt = 0.      #rho**2*t
    
    n = zeros(Int64, N)
    diff = zeros(Float64, N)
    dt = 0
    t_cluster = []
    dim = []
    
    while (t < T)
        
        #probabilities
        qD = rates.*map(length, sd)
        qI = kI*length(el)
        #qE = kE*length(over_thrE)
        q =  vcat(qD,qI)#,qE)
        qs = cumsum(q)
        qs /= qs[end]  #normalization

        #println()
        #println("Time: ", t)
        #println("-----------------------------------------------------------------------------")
        #println("Not normalized probability of extraction: \n", qE, "\n")
        #println("Not normalized probability of diffusion: \n", qD, "\n")
        #println("Not normalized probability of insertion: \n", qI, "\n")
        #println("All not normalized prob: diffusion, insertion, extraction \n", q ,"\n")
        #println("Normalized probabilities: \n", qs, "\n")

        # 1) it chooses the type of move
        move, type_of_move = choose_move(qs)  
        #println("type of move: ", type_of_move)
        dt = -log(rand(Float64))/sum(q)
        t = t + dt
        #=--------------------------------------DIFFUSION-------------------------------------------------=#
        if type_of_move == "diffusion"
            pc, rc, eta = diffusion(move, sd, eta)
            #println("The particle goes from ", pc, " to ", rc, "\n")
            rho_tot = rho_tot                          #diffusion does not change the density of the system
            #rho_i = rho_i

            lnm = unique(vcat(nn[pc,:],nn[rc,:]))                                # array of sites to update
            #println("lnm diff")
            m_loc = create_m_loc(lnm, eta)                #the same for insertion, diffusion and extraction
            #println("m_loc diff \n", m_loc,"\n")
            sd, isd = update_local_sd_isd(lnm, m, m_loc, sd, isd) # the same for insertion, diffusion and extraction
            el, pel = update_local_el_pel_diffusion(pc, rc, el, pel)
            lcc, cc, icc, over_thrE, over_soglia, dimeri, lccS, iccS, nₚ, nd, Nd, ns,tin, τᵤ, Ndis = update_local_lcc_cc_icc_diffusion(lnm, lcc, cc, icc, over_thrE, over_soglia, dimeri, lccS, iccS, eta, pc, rc, nₚ, nd, Nd, ns, tin, t, τᵤ, Ndis)
            #if (dimeri != Any[]) && (dimeri ∉ [dim[i][2] for i = 1:length(dim)])
            #    dim = push!(dim, [t, dimeri])
            #end
            
            m[lnm,:] = m_loc
        #=--------------------------------------INSERTION-------------------------------------------------=#
        elseif  type_of_move == "insertion"
            nc, col, eta = insertion(el, eta)  

            Nin += 1
            rho_tot +=  1. /N
            #rho_i[col] += 1. /N

            lnm = vcat(nc,nn[nc,:]) # array of sites to update
            #println("lnm ins")
            m_loc = create_m_loc(lnm, eta) #the same for insertion, diffusion and extraction
            #println("m_loc ins \n", m_loc,"\n")
            sd, isd =update_local_sd_isd(lnm, m, m_loc, sd, isd) # the same for insertion, diffusion and extraction
            el, pel = update_local_el_pel_insertion(nc, el, pel)
            lcc, cc, icc, over_thrE, over_soglia, dimeri, lccS, iccS, nₚ, nd, Nd , ns, tin = update_local_lcc_cc_icc_insertion(lcc, cc, icc, over_thrE, over_soglia, dimeri, lccS, iccS, eta, nc, nₚ, nd, Nd, ns, tin, t)

            #if (dimeri != Any[])  && (dimeri ∉ [dim[i][2] for i = 1:length(dim)])
            #    dim = push!(dim, [t, dimeri])
            #end
            
            m[lnm,:] = m_loc
        #=--------------------------------------EXTRACTION------------------------------------------------=#
        #= elseif  type_of_move == "extraction"
            ceot, col, eta = extraction(over_thrE, eta)
            
            over_soglia = push!(over_soglia, lcc[ceot])
            
            for i in lcc[ceot]
                n[i] += 1
            end
            n

            rho_tot -= length(lcc[ceot])/N    #change the density after an extraction
            #rho_i[col] -= length(lcc[ceot])/N

            # array of sites to update
            lnm = unique(nn[lcc[ceot],:]) 
            #println("lnm ext")
            m_loc = create_m_loc(lnm, eta) #the same for insertion, diffusion and extraction
            #println("m_loc ext \n", m_loc,"\n")
            sd, isd =update_local_sd_isd(lnm, m, m_loc, sd, isd) # the same for insertion, diffusion and extraction
            el, pel = update_local_el_pel_extraction(ceot, el, pel, lcc)
            lcc, cc, icc, over_thrE, lccS, iccS, nₚ, nd, Nd, ns, Next, next, tin, τₚ= update_local_lcc_cc_icc_extraction(
                                           ceot, lcc, cc, icc, over_thrE, lccS, iccS, nₚ, nd, Nd, ns, Next, next, tin, t, τₚ)
            m[lnm,:] = m_loc =#
        end
    
        densit = push!(densit, rho_tot)
        if ncol == 1
           cont, sum_dt, sum_rho_dt, sum_rho2_dt, Next, Ndis, Nin = statistics_1species(eta, t, rho_tot, nₚ, nd, Nd, ns, Next, next, cont, dt, sum_dt, sum_rho_dt, sum_rho2_dt, Ndis, Nin)
        end 
        
        τₚ = 0.
        τᵤ = 0.
        #println("sum_τₚ: ", sum_τₚ)
       
    end
    
    return densit, T
end

densit, T = tempo()
    

