function find_stable_time(sol::PowerGridSolution, powergrid::PowerGrid, fault, tol::Float64)
  n = length(sol.dqsol.u)
  pg = deepcopy(powergrid)
  pg = fault(pg)
  test_state = State(pg, sol.dqsol.u[1])
  idx_stop = -1
  for i in 1:n
    test_state = State(pg, sol.dqsol.u[i])
    # test_state.vec .= sol.dqsol.u[i]
    err = norm(err_from_root(test_state))
    # println("i=$i and err=$err")
    if err <= tol
      idx_stop = i
      break
    end
  end
  if idx_stop == -1
    return -Inf
  else
    t_stable = sol.dqsol.t[idx_stop]
    return t_stable
  end
end


function get_ok_lines(powergrid)
    ok = []
    mask = zeros(Bool, length(keys(powergrid.lines)))
    i = 0
    for lkey in keys(powergrid.lines)
        i += 1
        f = powergrid.lines[lkey].from
        t = powergrid.lines[lkey].to
        if typeof(powergrid.nodes[f]) == SwingEq && typeof(powergrid.nodes[t]) == SwingEq
            mask[i] = 0
            continue
        else
            push!(ok, lkey)
            mask[i] = 1
        end
    end
    return ok, mask
end