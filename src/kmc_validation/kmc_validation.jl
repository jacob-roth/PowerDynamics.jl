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

function load_kmc_dict(operatingdata_path::String, failure_idx::Int, num_seqs::Int=10) 
  Pnet = readdlm(operatingdata_path * "Pnet_$(failure_idx).txt")
  Qnet = readdlm(operatingdata_path * "Qnet_$(failure_idx).txt")
  xbar = readdlm(operatingdata_path * "xbar_$(failure_idx).txt")
  Vm = readdlm(operatingdata_path * "Vm_$(failure_idx).txt")
  Va = readdlm(operatingdata_path * "Va_$(failure_idx).txt")
  IDs = readdlm(operatingdata_path * "IDs.txt")
  ic_IDs = Int.(readdlm(operatingdata_path * "ic_IDs.txt"))

  Y_arr = Vector{Array{Complex}}(undef, num_seqs)
  for seq_idx in 1:num_seqs
    Y = readdlm(operatingdata_path * "Y_$(failure_idx)_seq_$(seq_idx).txt", '\t', ComplexF64)
    Y_arr[seq_idx] = Y
  end

  kmc_dict = Dict()
  kmc_dict[:Pnet] = Pnet
  kmc_dict[:Qnet] = Qnet
  kmc_dict[:xbar] = xbar
  kmc_dict[:Vm] = Vm
  kmc_dict[:Va] = Va
  kmc_dict[:Y] = Y_arr
  kmc_dict[:IDs] = IDs[:, failure_idx]
  kmc_dict[:ic_IDs] = ic_IDs

  return kmc_dict
end

function kmc2pd(fileout::String, opfmodeldata::Dict, kmc_dict::Dict,
                line_type::String, bus_type::String, seq_idx::Int,
                nobus_Bshunt::Bool=false, phys::Dict=physDefault)
  Y = kmc_dict[:Y][seq_idx]
  if isa(Y, AbstractArray{<:Complex})
    opfmodeldata[:Y] = imag.(Y)
  else
    opfmodeldata[:Y] = Y
  end

  optimal_values = Dict()
  optimal_values[:Vm] = kmc_dict[:Vm][seq_idx,:]
  optimal_values[:Va] = kmc_dict[:Va][seq_idx,:]
  optimal_values[:Pnet] = kmc_dict[:Pnet][seq_idx,:]
  optimal_values[:Qnet] = kmc_dict[:Qnet][seq_idx,:]
  optimal_values[:IDs] = kmc_dict[:IDs][seq_idx]
  optimal_values[:ic_IDs] = kmc_dict[:ic_IDs][seq_idx,:]

  return opf2pd(fileout, optimal_values, opfmodeldata, line_type, bus_type, nobus_Bshunt, phys)
end
