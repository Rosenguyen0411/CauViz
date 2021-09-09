#= Visualize the hierarchical structure given the position of nodes 
returned by function getpos =#

    function hierviz(Pos)
    
    Val = sum(convert(Array, value.(Pos)), dims = 1)

    #nNodes = Int(sqrt(length(Pos)))

    opt = round.((Val .+ 1 .+ nNodes) ./2)
    
    newNodeList = collect(1:nNodes)
    
    for i in 1:nNodes
        newNodeList[Int(opt[i])] = nodeList[i]
    end 
    
    
    Drawing(1800, 1000)
    origin()
    background("white")
    sethue("black")
    fontsize(20)
    
        radius = 30
        xloca = (opt.-1).*3 .*radius
    
    
    for i in 1:nNodes
        u = Point(xloca[i], 0) 
        text("u$i", u, halign=:center, valign=:middle)
        circle.(u, radius, :stroke)
    end
    
    for i in 1:nPlate
        for j in 1:(nNodes - length(plateList[i]) + 1)
            if sort(newNodeList[j:(j + length(plateList[i]) - 1)]) == sort(plateList[i])
                randomhue()
                box(Point(xloca[newNodeList[j]] - radius - 20*rand(), - radius - 20*rand()), 
                    Point(xloca[newNodeList[j + length(plateList[i]) - 1]] + radius + 20*rand(), radius + 20* rand()), :stroke)
            end
        end
    end
    
    finish()
    preview()
    

end