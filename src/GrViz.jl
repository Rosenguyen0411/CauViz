
function t_g(d::DAG)
    gs = topological_sort(d.a)
      vars = names(gs, 1)
      io = string("digraph $(d.name) {\n")

    for var in vars
      for (ind, entry) in enumerate(gs[:, var])
        if entry == 1
          io = string(io, "  $(var) -> $(vars[ind]);\n")
        end
      end
    end
    io = string(io, "}\n")
    return io
  end


function grViz(dt)
    @rput dt
    R"library(DiagrammeR)"
    R"grViz(dt)"
    end
