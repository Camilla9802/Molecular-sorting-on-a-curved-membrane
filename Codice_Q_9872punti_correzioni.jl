using SparseArrays, LinearAlgebra, Random, MAT, PyPlot, CSV, Tables, CavityTools
Random.seed!(1234)

NN = matread("neighbors_9872punti.mat")
NN["t"]
N_ = convert(Matrix{Int64}, NN["t"])

P = matread("coordinate_9872punti.mat")
p = P["p"]

b = CSV.read("baricentro_9872punti.csv", Tables.matrix, header=1)
nn = CSV.read("NN_9872punti.csv", Tables.matrix, header=1)

A = CSV.read("Aree_9872punti.csv", Tables.matrix, header=1)


nn_ = [nn[i,:] for i in 1:size(nn,1)]
for i = 1:length(nn_)
    j = 1
    while j <= length(nn_[i])
        if nn_[i][j]==0
            popat!(nn_[i], j)
            j -= 1
            nn_[i]
        else 
            j += 1
        end
    end
    nn_[i]
end
nn_

function tempo()
    #main
    nodes = 9872                 #number of nodes (a square of a number for the square lattice)                      
    ncol = 1                      #number of protein types (colors)
    T = 1e7                         #final time
    g = 10                       #aggregation constant (g >1 homologous neighbours reduce the prob )
    kI = 1e-6                     #insertion coefficient
    kE = 1e20                      #extraction coefficient 
    thrE = 50                    #extraction threshold
    kD = 1
    t = 0.0                        #initial time
    rho = 0.0
    
    v = zeros(Int64, nodes)
    for j = 1:nodes
        v[j] = count(i->(i != 0), nn[j, :])
    end
    
    eta = rand(1:ncol, nodes).*(rand(nodes).<rho)
    
    A_media = 0.0
    for i = 1:length(A)
        A_media += A[i]
    end
    A_media = A_media/length(A)
    
    for i = 1:length(A)
        A[i] = A[i]/A_media
    end
    
    eta_nn = copy(nn)
    for k = 1:length(eta)
        for i = 1:size(eta_nn,1)
            for j = 1:size(eta_nn,2)
                if eta[k] == 1 && eta_nn[i,j] == k
                    eta_nn[i,j] = 1
                elseif eta[k] == 0 && eta_nn[i,j] == k
                    eta_nn[i,j] = 0
                end
                eta_nn[i,j]
            end
        end
    end
    eta_nn
    
    function counting(eta::Vector{Int64})
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
        rho_tot = sum(nums[2:end])/nodes
    end
    #println("nums " ,nums)
    rho_i = zeros(Float64, 1, ncol)

    for (ind, col) in enumerate(nums[2:end])
        #println(ind,"  ", col)
        rho_i[ind] = nums[2:end][ind]/nodes
    end
    rho_tot, rho_i
    
    kD_eff = zeros(nodes)
    g_eff = zeros(nodes)
    n_tot = 6
    for i = 1:nodes
        kD_eff[i] = kD*n_tot/v[i]  
        g_eff[i] = g^(n_tot/v[i])  
    end
    
    eta_ = [eta[nn_[i]] for i = 1:length(nn_)]
    ne = (v - [sum(eta_[i].>0) for i=1:length(eta_)]).*(eta.>0)  #array of empty neighbors for each node
    eta1 = repeat(reshape(eta, nodes ,1), 1, 7)
    nh = sum((eta1 .== eta_nn) .& (eta1.>0) ,dims = 2)  #array of homologous neighbors for each node
    m = hcat(nh,ne)
    
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
                    #Each sub-array contains the position of the particles belonging to that cluster into i =#
    cc = -ones(Int64, nodes)     # it indicates the cluster (sub-array) where it is the particle
    icc = -ones(Int64, nodes)    # it gives us the position of the particle inside the sub-arrays 

    lccS = Vector{Int64}[] #perimeters of the connected components
    iccS = -ones(Int64, nodes) # it gives us the position of the particle inside the sub-arrays
    tin = Float64[]
    for i in 1:nodes
        if (eta[i] != 0) & (cc[i] == -1)
            #println(i)
            pos += 1
            cc[i] = pos
            icc[i] = 1
            push!(lcc, [i])
            l = length(lcc)
            if sum(eta_[i] .== eta[i]) != v[i]  #particle on the surface
                iccS[i] = 1
                push!(lccS, [i])
            else
                push!(lccS, [])
            end
            @inbounds for j in lcc[l]      #si prende un nodo del cluster l
                @inbounds for k in nn_[j]  #per il nodo j si prende un nodo vicino
                    if (eta[k] == eta[i]) & (cc[k] == -1)
                        #print(nn[j,:])
                        cc[k] = pos
                        push!(lcc[l], k)
                        icc[k] = length(lcc[l]) 

                        if sum(eta_[k] .== eta[i]) != v[i]   #particle on the surface
                        push!(lccS[l], k)
                        iccS[k] = length(lccS[l])
                        end 
                    end  
                end
            end
            if length(lcc[l]) > 1
                push!(tin, t)
                Nd += 1
                nd += length(lcc[l])/nodes  #num of particles inside a cluster
                ns += length(lccS[l])/nodes #number of particles in the surface of a cluster
            else
                nₚ += 1/nodes
                push!(tin, -1.)
            end
        end    
    end

    #clusters over size
    over_thrE = findall(x->length(x)>=thrE,lcc)
    
    function create_m_loc(lnm::Vector{Int64}, eta::Vector{Int64})
        ne_loc = [] #zeros(Int64, length(lnm))
        nh_loc = [] #zeros(Int64, length(lnm))
        for i in lnm
            push!(ne_loc, count(i->(i == 0), eta_[i]))
            push!(nh_loc, count(i->(i != 0), eta_[i]))
        end
        return m_loc = hcat(nh_loc, ne_loc)
    end
    
    function update_local_lcc_cc_icc_insertion(lcc::Vector{Vector{Int64}}, cc::Vector{Int64}, icc::Vector{Int64}, over_thrE::Vector{Int64}, lccS::Vector{Vector{Int64}}, iccS::Vector{Int64}, eta::Vector{Int64}, nc::Int64, nₚ::Float64, nd::Float64, Nd::Int64, ns::Float64, tin::Vector{Float64}, t::Float64, nn_=nn_::Array{Array{Int64, 1}, 1}, thrE::Int64=thrE, nodes::Int64=nodes)
        if eta[nc] ∉ eta_[nc] #1) isolated particle (monomer) 
            #println("homologous neighbors = 0")
            push!(lcc, [nc])
            push!(tin, -1)
            cc[nc] = length(lcc)
            icc[nc] = 1
            push!(lccS, [nc])
            iccS[nc] = 1
            nₚ += 1/nodes 

        elseif sum(eta[nc] .== eta_[nc]) == 1     # 2) homologous particle in the neighbourhood = 1
            #println("homologous neighbors = 1")
            ## I want to find the neighbour cluster
            plcc = cc[nn_[nc]][eta_[nc] .== eta[nc]][1]  #position in lcc
            ccn = copy(lcc[plcc])  #ccn = connected component neighbours
            #print("ccn: ", ccn)
            append!(lcc[plcc], nc)
            icc[nc] = length(lcc[plcc])
            cc[nc] = plcc
            if length(ccn) == 1 #join of two monomers
                #println("cluster: ", lcc[plcc])
                tin[plcc] = t    #the time where a new domain is formed
            end       

            nn1 = [nn[nc,i] for i = 1:length(nn[nc,:]) if nn[nc,i]!=0]
            an = nn1[eta_[nc] .== eta[nc]][1]
            #an = nn_[nc,eta_[nc] .== eta[nc]][1]                              #  
            #println("adjecient node to nc (an): ", an)                       #  
            if sum(eta_[an] .== eta[nc]) == v[nc] #(???)  #particle in the interior  #
                lccS[cc[an]][iccS[an]] = nc                                   #   
                iccS[nc] = iccS[an]                                           #
                iccS[an] = -1                                                 #
            else                                                              #    
                append!(lccS[plcc], nc)                                       #  perimeter 
                iccS[nc] = length(lccS[plcc])                                 # 
            end                                                               #  

            #update over_thrE
            if (length(lcc[plcc]) >= thrE) & (plcc ∉ over_thrE)
                push!(over_thrE, plcc)
            end
            #uptade nₚ, Nd and nd
            #print("ccn: ", ccn)
            #print("lcc[plcc]:",lcc[plcc])
            if length(ccn) == 1
                nₚ -= 1/nodes
                Nd += 1
                nd += length(lcc[plcc])/nodes
                ns += length(lccS[plcc])/nodes  #forse è più efficiente con 2/N
            else
                nd += 1/nodes
                if sum(eta_[an]) != v[nc] #(???)  
                    ns += 1/nodes                
                end
            end
            #print("np: ", nₚ, "   nd: ", nd, "  Nd: ", Nd)

        else
            #3) the new particle can marge two or more clusters (nh[nc] >= 2)
            #println("homologous neighbors => 2")
            cn = cc[nn_[nc][eta_[nc] .== eta[nc]]] #which clusters belonging to the homologous neighbors
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
            if sum(eta_[nc] .== eta[nc]) != v[nc] 
                push!(merged_clS, nc)
            end
            times_rm = Float64[]
            @inbounds for cl in clusters
                if length(lcc[cl]) >= thrE
                    over_thrE[cl] = over_thrE[end]
                    pop!(over_thrE)
                end
                if length(lcc[cl]) > 1      #update nd, nₚ 
                    nd -= length(lcc[cl])/nodes
                    ns -= length(lccS[cl])/nodes
                else
                    nₚ -= 1/nodes
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
            eta2 = [eta[nn_[nn_[nc]][i]] for i = 1:length(nn_[nn_[nc]])]
            p_in = nn_[nc] .* ([sum(eta2[i]) for i = 1:length(eta2)] .== v[nc])
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
            nd += length(merged_cl)/nodes
            ns += length(merged_clS)/nodes
            if num_dom == 0 
                Nd += 1
            elseif num_dom > 1
                Nd -= (num_dom -1) 
            end
            if length(merged_cl) >= thrE    
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
        return lcc, cc ,icc, over_thrE, lccS, iccS, nₚ, nd, Nd, ns, tin
    end
    
    function update_local_lcc_cc_icc_diffusion(lnm::Vector{Int64}, lcc::Vector{Vector{Int64}}, cc::Vector{Int64}, icc::Vector{Int64}, over_thrE::Vector{Int64}, lccS::Vector{Vector{Int64}}, iccS::Vector{Int64}, eta::Vector{Int64}, pc::Int64, rc::Int64, nₚ::Float64, nd::Float64, Nd::Int64, ns::Float64, tin::Vector{Float64}, t::Float64, τᵤ::Float64, Ndis::Int64, nn_=nn_::Array{Array{Int64, 1}, 1}, thrE::Int64=thrE)
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
                nd -= length(lcc[i])/nodes
                ns -= length(lccS[i])/nodes
                Nd -= 1
            else
                nₚ -= length(lcc[i])/nodes
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
        iccS[nn_[rc]] .= -1 
        #println("nₚ: ", nₚ)
        #println("nd: ", nd)
        #println("Nd: ", Nd)


        #add new clusters
        tlᵤ = Float64[]
        colors = Int64[]
        starting_points = [k for k in nn_[pc] if eta[rc] == eta[k]]
        @inbounds for s in starting_points
            if cc[s] ∉ colors
                c = length(lcc) + 1
                cc[s] = c
                icc[s] = 1
                push!(colors, c)
                push!(lcc, [s])          
                if sum(eta_[s] .== eta[s]) != v[s]  #particle on the surface            
                    push!(lccS, [s])
                    iccS[s] = 1
                end
                @inbounds for j in lcc[c]
                    @inbounds for k in nn_[j]
                        if (eta[k] == eta[j]) & (cc[k] ∉ colors)
                            push!(lcc[c], k)
                            if sum(eta_[k] .== eta[k]) != v[k]
                                push!(lccS[c], k)
                                iccS[k] =length(lccS[c])
                            end
                            cc[k] = c
                            icc[k] =length(lcc[c]) 
                        end
                    end               
                end
                if (length(lcc[c]) >= thrE) & (length(lcc) ∉ over_thrE)         
                    push!(over_thrE, length(lcc))
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
                    nd += length(lcc[c])/nodes  #num of particles inside domains
                    ns += length(lccS[c])/nodes            
                    #println("Nd: ", Nd)
                else
                    nₚ += 1/nodes
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

        return lcc, cc, icc, over_thrE, lccS, iccS, nₚ, nd, Nd, ns, tin, τᵤ, Ndis            
    end
    
    function update_local_lcc_cc_icc_extraction(ceot::Int64, lcc::Vector{Vector{Int64}}, cc::Vector{Int64}, icc::Vector{Int64}, over_thrE::Vector{Int64}, lccS::Vector{Vector{Int64}}, iccS::Vector{Int64}, nₚ::Float64, nd::Float64, Nd::Int64, ns::Float64, Next::Int64, next::Float64, tin::Vector{Float64}, t::Float64, τₚ::Float64)
    
        nd -= length(lcc[ceot])/nodes #update nd
        ns -= length(lccS[ceot])/nodes
        Nd -= 1
        Next += 1 
        next += length(lcc[ceot])/nodes

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
    end
    
    function statistics_1species(eta,t, rho_tot, nₚ, nd, Nd, ns, Next, next, cont, dt, sum_dt, sum_rho_dt, 
                                 sum_rho2_dt, Ndis, Nin, T=T, nodes=nodes, thrE=thrE, g=g, ncol=ncol, kI=kI; points = 600000)

        sum_dt += dt
        sum_rho_dt += dt*rho_tot
        sum_rho2_dt += dt*rho_tot^2
        cont += 1

        if cont == points
            avg_rho = sum_rho_dt/sum_dt
            avg_rho2 = sum_rho2_dt/sum_dt
            var = avg_rho2 - avg_rho^2
            
            write_on_file(eta, t, avg_rho, var, Next, Nin, Ndis, n, densit)
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
    
    function write_on_file(eta, t, avg_rho, var, Next, Nin, Ndis, n, densit, T::Float64=T, nodes::Int64=nodes, thrE::Int64=thrE, g::Int64=g, ncol::Int64=ncol, kI::Float64=kI)
        #data series
        mol_sorting = open("Mol_sorting3(nodes_$(nodes)_g_$(g)_kI_$(kI)_thrE_$(thrE))", "a")
        data = "$t" *"    " * "$avg_rho" * "    " * "$var" * "    " * "$Next" * "    " * "$Nin" * "    " * "$Nd" * "\n"
        write(mol_sorting, data)
        close(mol_sorting)
        
        system = open("Extraction3(nodes_$(nodes)_g_$(g)_kI_$(kI)_thrE_$(thrE))", "w")
        write(system, "$n \n")
        close(system)
        
        corr = open("Correlation3(nodes_$(nodes)_g_$(g)_kI_$(kI)_thrE_$(thrE))", "w")
        data1 = "$t_cluster" * "\n"
        write(corr, data1)
        close(corr)
        
        #diffusion = open("Diffusion_corr(nodes_$(nodes)_g_$(g)_kI_$(kI)_thrE_$(thrE))", "w")
        #write(diffusion, "$diff \n")
        #close(diffusion)
        
        #conf = open("State(nodes_$(nodes)_g_$(g)_kI_$(kI)_thrE_$(thrE))", "a")
        #write(conf, "$eta \n")
        #close(conf)
        
        #density = open("Density(nodes_$(nodes)_g_$(g)_kI_$(kI)_thrE_$(thrE))", "a")
        #write(density, "$avg_rho \n")
        #close(density)
        
    end
    
    #statistics
    cont = 0
    sum_dt = 0.          
    sum_rho_dt = 0.       #rho*t
    sum_rho2_dt = 0.      #rho**2*t
    
    n = zeros(Int64, nodes)
    diff = zeros(Float64, nodes)
    dt = 0
    densit = []
    t_cluster = []
    
    Q = ExponentialQueue(2*nodes + 1)
    for j = 1:nodes
        if eta[j] != 0
            Q[j] = kD * m[j,2]/g^m[j,1]
        else
            Q[j + nodes] = kI * A[j]
        end
    end

    Q[2*nodes + 1] = kE*length(over_thrE)
    
    while (t < T)
        
        ind, dt = pop!(Q)
        #println("Time: $t \n")
        #println("i: $i" * "   " * "dt: $dt" * "\n")
        if ind <= nodes   #DIFFUSION
            pc = ind                            # 2) random choice of the particle to move
            rc = rand(nn_[pc][findall(x->x==0, eta_[pc])]) # 3) select an empty neighbors of pc at random
            
            #println("The particle goes from ", pc, " to ", rc, "\n")
            
            prob = v[pc]/v[rc]
            if rand() < prob
                #update the system
                eta[rc] = eta[pc]       #move particle pc in site rc
                eta[pc] = 0
                eta_[pc]=eta[nn_[pc]]
                eta_[rc]=eta[nn_[rc]]
                for neig in nn_[pc]
                    eta_[neig]=eta[nn_[neig]]
                end
                for neig in nn_[rc]
                    eta_[neig]=eta[nn_[neig]]
                end
                
                for i = 1:length(eta)
                    if eta[i] == 1
                        diff[i] += dt
                    end
                end
                diff

                lnm = unique(vcat(nn_[pc],nn_[rc]))                                # array of sites to update
                m_loc = create_m_loc(lnm, eta)                #the same for insertion, diffusion and extraction
                lcc, cc, icc, over_thrE, lccS, iccS, nₚ, nd, Nd, ns,tin, τᵤ, Ndis = update_local_lcc_cc_icc_diffusion(lnm, lcc, cc, 
                                                     icc, over_thrE, lccS, iccS, eta, pc, rc, nₚ, nd, Nd, ns, tin, t, τᵤ, Ndis)

                Q[2*nodes + 1] = kE*length(over_thrE)

                m[lnm,:] = m_loc
                for j in lnm
                    if eta[j] != 0
                        Q[j] = kD * m[j,2]/g^m[j,1]
                        Q[j + nodes] = 0
                    else
                        Q[j + nodes] = kI * A[j]
                    end
                end
                
            else
                
                lnm = unique(vcat(nn_[pc],nn_[rc]))
                for j in lnm
                    if eta[j] != 0
                        Q[j] = kD_eff[j] * m[j,2]/g_eff[j]^m[j,1]
                        Q[j + nodes] = 0
                    else
                        Q[j + nodes] = kI * A[j]
                    end
                end
            end
            
        elseif nodes < ind <= 2*nodes  #INSERTION
            nc = ind-nodes          # 2) choose a node at random (nc = chosen node)
            
            #println("The particle is inserted in $nc \n")
            
            col = rand(1:ncol)      # 3) choose the protein type to insert
            #update the system
            eta[nc] = col
            eta_[nc]=eta[nn_[nc]]
            for neig in nn_[nc]
                eta_[neig]=eta[nn_[neig]]
            end
            
            for i = 1:length(eta)
                if eta[i] == 1
                    diff[i] += dt
                end
            end
            diff
        
        
            Nin += 1
            rho_tot +=  1. /nodes
            #rho_i[col] += 1. /nodes
         
            lnm = vcat(nc,nn_[nc]) # array of sites to update
            m_loc = create_m_loc(lnm, eta) #the same for insertion, diffusion and extraction
            lcc, cc, icc, over_thrE, lccS, iccS, nₚ, nd, Nd , ns, tin = update_local_lcc_cc_icc_insertion(lcc, cc, icc, 
                                                            over_thrE, lccS, iccS, eta, nc, nₚ, nd, Nd, ns, tin, t)
            
            Q[2*nodes + 1] = kE*length(over_thrE)
        
            m[lnm,:] = m_loc
            
            for j in lnm
                if eta[j] != 0
                    Q[j] = kD * m[j,2]/g^m[j,1]
                end
            end
            
        else    #EXTRACTION
            ceot = rand(over_thrE)         #ceot = select a cluster from over_thrE (position of the cluster in sd)
            col = eta[lcc[ceot]][1]        #color of the cluster that has to be extracted 
            cluster = lcc[ceot]
            t_in = tin[ceot]
            
            #println("The cluster extracted in is $ceot \n")
            #update the system
            eta[lcc[ceot]] .= 0  #extraction
            for elemento in lcc[ceot]
                eta_[elemento]=eta[nn_[elemento]]
                for neig in nn_[elemento]
                    eta_[neig]=eta[nn_[neig]]
                end 
            end
            
            for i = 1:length(eta)
                if eta[i] == 1
                     diff[i] += dt
                end
            end
            diff
        
            for i in lcc[ceot]
                n[i] += 1
            end
            n
        
            rho_tot -= length(lcc[ceot])/nodes    #change the density after an extraction
            #rho_i[col] -= length(lcc[ceot])/nodes
        
            # array of sites to update
            lnm = unique(nn[lcc[ceot], :]); lnm = [lnm[i] for i = 1:length(lnm) if lnm[i]!=0] 
            m_loc = create_m_loc(lnm, eta) #the same for insertion, diffusion and extraction
            lcc, cc, icc, over_thrE, lccS, iccS, nₚ, nd, Nd, ns, Next, next, tin, τₚ= update_local_lcc_cc_icc_extraction(
                                           ceot, lcc, cc, icc, over_thrE, lccS, iccS, nₚ, nd, Nd, ns, Next, next, tin, t, τₚ)
            
            t_cluster = push!(t_cluster, [t_in, τₚ, cluster])
            
            Q[2*nodes + 1] = kE*length(over_thrE)
            
            m[lnm,:] = m_loc
            
            for j in lnm
                if eta[j] != 0
                    Q[j] = kD * m[j,2]/g^m[j,1]
                else
                    Q[j] = 0
                    Q[j + nodes] = kI * A[j]
                end
            end
            
        end
    
        densit = push!(densit, rho_tot)
        
        if ncol == 1
            cont, sum_dt, sum_rho_dt, sum_rho2_dt, Next, Ndis, Nin = statistics_1species(eta, t, rho_tot, nₚ, nd, Nd, ns, Next, next, cont, dt, sum_dt, sum_rho_dt, sum_rho2_dt, Ndis, Nin)
        end  
        
        τₚ = 0.
        τᵤ = 0.
        #println("sum_τₚ: ", sum_τₚ)
        #dt = -log(rand(Float64))/sum(q)
        t = t + dt
    end
    
    return densit, T
end

densit, T = tempo()