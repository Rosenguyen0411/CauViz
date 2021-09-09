## Get the optimal position of nodes given the list of nodes and plates

function getpos(N, plateList)

    nPlate = length(plateList)

    nNodes = 0
    for i in 1:length(N)
        nNodes += length(N[i])
    end

    nodeList = collect(1:nNodes)

    plateLevelList = []

    for i in 1: nPlate 
        h = []
        for k in 1:length(plateList[i])
            for j in 1: length(N)
                if plateList[i][k] ∈ N[j]
                    push!(h,j)
                end
            end
        end
        push!(plateLevelList, h)
    end

    rangeList=[]
    for i in 1: nPlate 
        h = []
        push!(h,minimum(plateLevelList[i]), maximum(plateLevelList[i]))
        push!(rangeList, h)
    end

    #model = Model(with_optimizer(Mosek.Optimizer, QUIET=false, INTPNT_CO_TOL_DFEAS=1e-7))
    model = Model(GLPK.Optimizer)

    nX = nNodes^2 - nNodes

    nA = 0
    for j in 1:nPlate
        for i in rangeList[j][1] : rangeList[j][2]
            nA += length(setdiff(N[i], plateList[j]))
        end
    end


    @variable(model, A[1:nA] >= 0)
    @variable(model, B[1:nA], Bin)
    @variable(model, X[1:nNodes,1:nNodes], Int)
    #@variable(model, Y[1:14,1:2], Bin)

    @objective(model, Min, sum(A[i] for i in 1:nA))


    w = 1
    for j in 1:nPlate
        for i in rangeList[j][1] : rangeList[j][2]
            for k in setdiff(N[i], plateList[j])
                @constraint(model, sum(X[k,plateList[j][h]] for h in 1:length(plateList[j])) <= 
                    A[w] + length(plateList[j]))
                @constraint(model, -sum(X[k,plateList[j][h]] for h in 1:length(plateList[j])) <= 
                    A[w] + length(plateList[j]))
                @constraint(model, sum(X[k,plateList[j][h]] for h in 1:length(plateList[j])) >=  
                    length(plateList[j]) - A[w] - 1000*B[w])
                @constraint(model, sum(X[k,plateList[j][h]] for h in 1:length(plateList[j])) <=  
                    1000 - length(plateList[j]) + A[w] - 1000*B[w])
                w = w + 1
            end
        end 
    end 

    for i in 1:(nNodes-1)
        for j in (i+1):nNodes
            @constraint(model, X[i,j] + X[j,i] == 0)
        end
    end



    for i in 1:nNodes
            @constraint(model, X[i,i] == 0)
    end


    for i in 1:nNodes
        for j ∈ nodeList[setdiff(1:nNodes, i)]
            for k ∈ nodeList[setdiff(1:nNodes, [i,j])]
                @constraint(model, X[i,j] + X[j,k] - X[i,k]  <= 1)
            end
        end
    end

    @constraint(model, -1 .<= X .<= 1)


    optimize!(model)
    @show value.(X)
    return(nNodes)
end
