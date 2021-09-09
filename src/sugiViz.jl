"""
    Assigns layers using the longest path method.
    Arguments:
    adj_list        Directed graph in adjacency list format
"""

function layerAssmtLongestpathRec(adj_list, layers, j)
    # Look for all children of vertex i, try to bump their layer
    n = length(adj_list)
    for i in adj_list[j]
        if layers[i] == -1 || layers[i] <= layers[j] #if the layer of immediate precedor node == -1 or lower than current node, add 1 to precedor layer
            layers[i] = layers[j] + 1
            layerAssmtLongestpathRec(adj_list, layers, i) # recursive on precedor node
        end
    end
end


function layerAssmtLongestpath(adj_list)
    n = length(adj_list)
    layers = fill(-1, n) #assign each node to a layer. each with value -1

    for j in 1:n #
        out_deg = 0
        for i in 1:n
            if j in adj_list[i]
                out_deg += 1
            end
        end
        if out_deg == 0 # node without out-edge will be assigned to layer 1
            # Start recursive walk from this vertex
            layers[j] = 1
            layerAssmtLongestpathRec(adj_list, layers, j)
        end
    end

    return layers
end





"""
    Given a layer assignment, introduce dummy vertices to break up
    long edges (more than one layer)
    Arguments:
    orig_adj_list   Original directed graph in adjacency list format
    layers          Assignment of original vertices
"""
function layerAssmtDummy(orig_adj_list, layers)
    adj_list = deepcopy(orig_adj_list)

    # This is essentially
    # for i in 1:length(adj_list)
    # but the adj_list is growing in the loop
    i = 1
    while i <= length(adj_list)
        for (k,j) in enumerate(adj_list[i]) # k is the counter and j is the value of adj_list[i] at kth position
            if layers[j] - layers[i] > 1 # if the difference between the layer is larger than 1 => long edge
                # Need to add a dummy vertex
                new_v = length(adj_list) + 1 
                adj_list[i][k] = new_v  # Replace the node of the long-edge with the dummy node
                push!(adj_list, Int[j])  # Make the node of the long-edge precedor of the dummy node
                push!(layers, layers[i]+1)  # Layer for the dummy node 
            end
        end
        i += 1
    end

    return adj_list, layers
end


"""
    Given a layer assignment, decide a permutation for each layer
    that attempts to minimizes edge crossings using the barycenter
    method proposed in the Sugiyama paper.
    Arguments:
    adj_list        Directed graph in adjacency list format
    layers          Assignment of vertices
    layer_verts     Dictionary of layer => vertices
"""
function orderingBarycentric(adj_list, layers, layer_verts)
    num_layers = maximum(layers)
    n = length(adj_list)

    for iter in 1:5
        # DOWN
        for L in 1:num_layers-1
            # Calculate barycenter for every vertex in next layer
            cur_layer = layer_verts[L]
            next_layer = layer_verts[L+1]
            barys = zeros(n)
            out_deg = zeros(n)
            
            for (p,i) in enumerate(cur_layer)
                # Because of the dummy vertices we know that
                # all vertices in adj list for i are in next layer
                for (q,j) in enumerate(adj_list[i])
                    barys[j] = barys[j] .+ p
                    out_deg[j] = out_deg[j] .+ 1
                end
            end
            
            barys = barys ./ out_deg
            
            # Arrange next layer by barys, in ascending order
            next_layer_barys = [barys[j] for j in next_layer]
           
            layer_verts[L+1] = next_layer[sortperm(next_layer_barys)]
        end
        
        # UP
        for L in num_layers-1:-1:1
            # Calculate barycenters for every vertex in cur layer
            cur_layer = layer_verts[L]
            next_layer = layer_verts[L+1]
            barys = zeros(n)
            out_deg = zeros(n)
            for (p,i) in enumerate(cur_layer)
                # Because of the dummy vertices we know that
                # all vertices in adj list for i are in next layer
                # We need to know their positions in next layer
                # though, unfortunately. Probably a smarter way
                # to do this step
                for (q,j) in enumerate(adj_list[i])
                    # Find position in next layer
                    for (r,k) in enumerate(next_layer)
                        if k == j
                            barys[i] = barys[i] .+ r
                            out_deg[i] = out_deg[i] .+ 1
                            break
                        end
                    end
                end
            end
            barys = barys ./ out_deg
            # Arrange cur layer by barys, in ascending order
            cur_layer_barys = [barys[j] for j in cur_layer]
            #println("UP $L")
            #println(cur_layer)
            #println(cur_layer_barys)
            layer_verts[L] = cur_layer[sortperm(cur_layer_barys)]
        end
        # Do something with phase 2 here - don't really understand
    end

    return layer_verts
end


"""
    Given a layer assignment, decide a permutation for each layer
    that minimizes edge crossings using integer programming.
    Based on the IP described in
      M. Junger, E. K. Lee, P. Mutzel, and T. Odenthal.
      A polyhedral approach to the multi-layer crossing minimization problem.
      In G. Di Battista, editor, Graph Drawing: 5th International Symposium,
      GD ’97, volume 1353 of Lecture Notes in Computer Science, pages 13–24,
      Rome, Italy, September 1997. Springer-Verlag.
    Arguments:
    adj_list        Directed graph in adjacency list format
    layers          Assignment of vertices
    layer_verts     Dictionary of layer => vertices (initial perm.)
    Returns:
    new_layer_verts An improved dictionary of layer => vertices (opt. perm.)
"""
function orderingIp(adj_list, layers, layer_verts)
    num_layers = maximum(layers)

    m = Model(GLPK.Optimizer)

    # Define crossing binary variables
    @variable(m, c[L=1:num_layers,        # for each layer
                 i=layer_verts[L],      # for each vertex in this layer
                 j=adj_list[i],         # and vertex in the next layer
                 k=layer_verts[L],      # for each vertex in this layer
                 l=adj_list[k]], Bin)   # and vertex in the next layer

    # Objective: minimize crossings
    @objective(m, Min, sum(c))

    # Define permutation variables for each layer
    # We'll define for both (i,j) and (j,i), and ensure they consistency
    # by adding constraints. Presolve in the IP solver will simplify
    # the problem for us by removing one of the variables.
    @variable(m, x[L=1:num_layers, i=layer_verts[L], j=layer_verts[L]], Bin)
    for L in 1:num_layers
        for i in layer_verts[L]
            for j in layer_verts[L]
                j <= i && continue  # Don't double-add
                # Ensure x[i,j] and x[j,i] are consistent
                @constraint(m, x[L,i,j] == 1 - x[L,j,i])
                # And ensure that triples are consistent
                for k in layer_verts[L]
                    k <= j && continue
                    @constraint(m, 0 <= x[L,i,j] + x[L,j,k] - x[L,i,k] <= 1)
                end
            end
        end
    end

    # Link permutations to crossings
    for L in 1:num_layers-1
        # For all (i,j)
        for i in layer_verts[L]
        for j in adj_list[i]
            # For all (k,l)
            for k in layer_verts[L]
            k == i && continue  # Can't cross if starting from same vertex!
            for l in adj_list[k]
                @constraint(m, -c[L,i,j,k,l] <= x[L+1,j,l] - x[L,i,k])
                @constraint(m,  c[L,i,j,k,l] >= x[L+1,j,l] - x[L,i,k])
            end
            end
        end
        end
    end

    # Solve the IP
    optimize!(m)

    # Extract permutation from solution
    x_sol = getvalue(x)
    new_layer_verts = [L => Int[] for L in 1:num_layers]
    for L in 1:num_layers
        old_perm = layer_verts[L]
        # For each vertex, count the number of times it is "in front"
        # of the other vertices. The higher this number, the earlier
        # the vertex appears in the layer.
        scores = zeros(length(old_perm))
        for (p,i) in enumerate(old_perm)
            for j in old_perm
                i == j && continue
                if round(Integer, x_sol[L,i,j]) == 1
                    # i appears before j
                    scores[p] += 1
                end
            end
        end
        new_layer_verts[L] = old_perm[sortperm(scores,rev=true)]
    end

    return new_layer_verts
end


########################################################################


"""
    Given a layer assignment and permutation, decide the coordinates for
    each vertex. The objective is to encourage straight edges, especially
    for longer edges. This function uses an integer program to decide the
    coordinates (although it is solved as a linear program), as described in
      Gansner, Emden R., et al.
      A technique for drawing directed graphs.
      Software Engineering, IEEE Transactions on 19.3 (1993): 214-230.
    Arguments:
    adj_list        Directed graph in adjacency list format
    layers          Assignment of vertices
    layer_verts     Dictionary of layer => vertices (final perm.)
    orig_n          Number of original (non-dummy) vertices
    widths          Width of each vertex
    xsep            Minimum seperation between each vertex
    Returns:
    layer_coords    For each layer and vertex, the x-coord
"""
function coordIp(adj_list, layers, layer_verts, orig_n, widths, xsep)
    num_layers = maximum(layers)

    m = Model(GLPK.Optimizer)
    

    # One variable for each vertex
    @variable(m, x[L=1:num_layers, i=layer_verts[L]] >= 0)

    # Constraint: must respect permutation, and spacing constraint
    for L in 1:num_layers
        for i in 1:length(layer_verts[L])-1
            a = layer_verts[L][i]
            b = layer_verts[L][i+1]
            @constraint(m, x[L,b] - x[L,a] >=
                (widths[a] + widths[b])/2 + xsep) #the distance between two
                #adjacent node must be at least xsep
        end
    end

    # Objective: minimize total misalignment
    # Use the weights from the Ganser paper:
    #   1 if both nodes "real"
    #   2 if one of the nodes is "real"
    #   8 if neither node is "real"
    # We use absolute distance in the objective so we'll need
    # auxilary variables for each pair of edges
    obj = AffExpr()
    
    @variable(m, absdiff[L=1:num_layers-1,
                        i=layer_verts[L], j=adj_list[i]] >= 0)
    
    for L in 1:num_layers-1
        for i in layer_verts[L]
            for j in adj_list[i]
                #add 2 constraints to deal with absolute value in objective
                @constraint(m, absdiff[L,i,j] >= x[L,i] - x[L+1,j]) 
                @constraint(m, absdiff[L,i,j] >= x[L+1,j] - x[L,i])
                if i > orig_n && j > orig_n
                    # Both dummy vertices
                    obj += 8*absdiff[L,i,j]
                elseif (i <= orig_n && j >  orig_n) ||
                       (i >  orig_n && j <= orig_n)
                    # Only one dummy vertix
                    obj += 2*absdiff[L,i,j]
                else
                    # Both real
                    obj += absdiff[L,i,j]
                end
            end
        end
    end
    @objective(m, Min, obj)

    # Solve it...
    optimize!(m)

    # ... and mangle the solution into shape
    x_sol = getvalue.(x)
    locs_x = zeros(length(layers))
    for L in 1:num_layers
        for i in layer_verts[L]
            locs_x[i] = x_sol[L,i]
        end
    end
    return locs_x
end


########################################################################
function layoutTree(adj_list,
    labels::Vector;
    filename    = "",

    cycles      = true,
    ordering    = :optimal,
    coord       = :optimal,
    xsep        = 50,
    x_start       = 100,
    ysep        = 100,
    scale       = 0.05,
    labelpad    = 1.2,
    background  = "white")
# Calculate the original number of vertices
n = length(adj_list)

# 1     Cycle removal
if cycles
# Need to remove cycles first
error("Cycle removal not implemented!")
end

# 2     Layering
# 2.1   Assign a layer to each vertex
layers = layerAssmtLongestpath(adj_list)
num_layers = maximum(layers)
# 2.2  Create dummy vertices for long edges
adj_list, layers = layerAssmtDummy(adj_list, layers)
orig_n, n = n, length(adj_list)


# 3     Vertex ordering [to reduce crossings]
# 3.1   Build initial permutation vectors
layer_verts = Dict([L => Int[] for L in 1:num_layers])

for i in 1:n
push!(layer_verts[layers[i]], i)
end


# 3.2  Reorder permutations to reduce crossings
if ordering == :barycentric
    layer_verts = orderingBarycentric(adj_list, layers, layer_verts)
elseif ordering == :optimal
    layer_verts = orderingIp(adj_list, layers, layer_verts)
end


# 4     Vertex coordinates [to straighten edges]
# 4.1   Place y coordinates in layers
locs_y = zeros(n)

# Luxor origin (0,0) is the upper-left corner. Lower layer has higher y location.
rev_y = sort(1:num_layers, rev = true) #reverse the layer
#ysep = 100

for L in 1:num_layers
    for (x,v) in enumerate(layer_verts[L])
        locs_y[v] = (rev_y[L])*ysep
    end
end

# 4.2   Get widths of each label, if there are any
widths  = ones(n)
widths[orig_n+1:n]  .= 0 #dummy node has 0 width
heights = ones(n)
heights[orig_n+1:n] .= 0 #dummy node has 0 height

# Note that we will convert these sizes into "absolute" units
# and then work in these same units throughout. The font size used
# here is just arbitrary, and unchanging. This hack arises because it
# is meaningless to ask for the size of the font in "relative" units
# but we don't want to collapse to absolute units until the end.

if length(labels) == orig_n
    extents = textextents.(labels)
    for (i,(width,height)) in enumerate(extents)
        widths[i]  = reduce(hcat, extents)[3,i]
        heights[i] = reduce(hcat, extents)[4,i]
    end
end



locs_x = coordIp(adj_list, layers, layer_verts, orig_n, widths, xsep) .+ x_start #= add x_start to shift all x location
to the right by 100 =#

# 4.3   Summarize vertex info
max_x, max_y = maximum(locs_x), maximum(locs_y)
max_w, max_h = maximum(widths), maximum(heights)

## find x_pad to the right and add the pad to the border of the most right node 
x_pad = minimum(locs_x - widths/2)
w = ceil.(max_x .+ widths[locs_x .== max_x]/2 .+ x_pad)[1] # width of the drawing is the max_x + width of the node + x_pad
h = ceil(max_y + max_h/2 + ysep) # height of the drawing is the max_y + max_height + y_pad

# 5     Draw the tree

#=
@draw begin
    Luxor.fontsize(20)
    fontface("Noto-Sans") 

    #widths = reduce(hcat, textextents.(labels))[3,:]
    widths = ceil.(widths .* 0.6)
    #push!(widths, 0)

    #heights = reduce(hcat, textextents.(labels))[4,:]
    #push!(heights, 0)
    heights = repeat([max_h],n)
    heights = ceil.(heights)
    

    for i in 1:orig_n
        #origin()
        point = Luxor.Point(locs_x[i], locs_y[i])
        
        squircle(point, widths[i], heights[i], rt = 0.3, :stroke)
        #textcentered(labels[i], point)
        Luxor.text(labels[i], point, halign=:center, valign=:middle)
    end 

    for L in 1:num_layers, i in layer_verts[L], j in adj_list[i]
        if i <= orig_n #if the destination node is a dummy then dont add arrow
            Luxor.arrow(Point(locs_x[j], locs_y[j] + heights[j]), Point(locs_x[i], locs_y[i] - heights[i]), linewidth=2)
        else 
            Luxor.line(Point(locs_x[j], locs_y[j] + heights[j]), Point(locs_x[i], locs_y[i] + heights[i]), :stroke)
        end
    end

    end
    =#

Drawing(w, h, "test.svg")

Luxor.background("$background")
#sethue("red")
Luxor.fontsize(20)
fontface("Noto-Sans")

for i in 1:orig_n
    
    point = Luxor.Point(locs_x[i], locs_y[i])
    
    squircle(point, widths[i]/1.5, ceil(max_h), rt = 0.3, :stroke)
    
    Luxor.text(labels[i], point, halign=:center, valign=:middle)
end 

for L in 1:num_layers, i in layer_verts[L], j in adj_list[i]
    if i <= orig_n #if the destination node is a dummy then dont add arrow
        Luxor.arrow(Point(locs_x[j], locs_y[j] + ceil(max_h)), Point(locs_x[i], locs_y[i] - ceil(max_h)), linewidth=2)
    else 
        Luxor.line(Point(locs_x[j], locs_y[j] + ceil(max_h)), Point(locs_x[i], locs_y[i] + ceil(max_h)), :stroke)
    end
end

finish()
preview()

end

