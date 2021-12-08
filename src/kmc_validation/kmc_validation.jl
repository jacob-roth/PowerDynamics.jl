function validate(op_kmc::State,
                  op_l_kmc::State,
                  degraded_powergrid::PowerGrid,
                  # ls::Array{I,1},
                  l::I,
                  tol::Float64=1e-5
                  ) where I <: Integer
  tfault_start = 1.0
  tfault_end = Inf
  tsim_start = 0.0
  tsim_end = Inf
  tspan_fault = (tfault_start, tfault_end)
  tspan_sim = (tsim_start, tsim_end)
  # multi_line_fault = LineFailures(line_names=["branch$(l)" for ll in ls], tspan_fault=tspan_fault);
  
  ##
  ## stability in toplevel degraded network
  ##   --> verify that `op_kmc`` is stable and compare to `op_PD`
  ##
  op_PD = find_operationpoint(degraded_powergrid, sol_method=:dynamic)
  nofault = ChangeInitialConditions(node = "bus1", var = :ω, f = Inc(0.0)); # do nothing
  solution_op_kmc = simulate(nofault, degraded_powergrid, op_kmc, timespan); # simulate PD from `op_kmc`
  tstable_op_kmc = find_stable_time(solution_op_kmc, tol)
  
  ##
  ## stability in toplevel degraded network + line `l` failure
  ##   --> verify that `op_l_kmc`` is stable and compare to `op_l_PD`
  ##
  single_line_fault = LineFailure(line_name="branch$(l)", tspan_fault=tspan_fault); # fail line `l`
  solution_op_l_PD = simulate(nofault, degraded_powergrid, op_kmc, timespan); # simulate PD from `op_kmc` to obtain `op_l_PD`
  tstable_op_kmc = find_stable_time(solution_op_l_PD, tol)
  op_l_diff = norm(solution_op_l_PD.dqsol.u[end] .- op_l_kmc.vec)
  
end

function find_stable_time(sol::PowerGridSolution, tol::Float64)
  (:u, :u_analytic, :errors, :t, :k, :prob, :alg, :interp, :dense, :tslocation, :destats, :retcode)
  n = length(sol.dqsol.u)
  test_state = State(sol.powergrid, sol.dqsol.u[1])
  idx_stop = -1
  for i in 1:n
    test_state = State(sol.powergrid, sol.dqsol.u[i])
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

function kmc2pd(M::Dict)
  """
  tbd
  """
  Vm_internal = M[:Vm_internal]
  Va_intenal = M[:Va_intenal]
  Vm_IDs = M[:Vm_IDs]
  Va_IDs = M[:Va_IDs]
  Vm = M[:Vm]
  Va = M[:Va]
  Y = M[:Y]
  Pnet = M[:Pnet]
  Qnet = M[:Qnet]
  
  power_injections = M[:power_injections]    
  Pg = M[:Pg]
  Qg = M[:Qg]
  Pd = M[:Pd]
  Qd = M[:Qd]
  save_xbar = M[:save_xbar]
end

function load_kmc_output(filename::String)
  """
  tbd
  """
    
filename = 
end