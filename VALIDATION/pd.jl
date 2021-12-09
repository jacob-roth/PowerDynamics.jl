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