using PowerDynamics: read_powergrid, Json, Inc, find_operationpoint, ChangeInitialConditions, LineFailure, PowerPerturbation, simulate

##
## 118bus with swingeq
##

powergrid = read_powergrid(joinpath(@__DIR__, "118notap.json"), Json)
operationpoint = find_operationpoint(powergrid, sol_method = :dynamic) # solve_powerflow=true, 

operatingdata_path = "/Users/jakeroth/git/PowerDynamics.jl/examples/118bus/"
imported_vec = import_own_operatingpoint(powergrid, operatingdata_path)
imported_operationpoint = State(powergrid, imported_vec)
norm(imported_operationpoint.vec - operationpoint.vec)

timespan = (0.0, 25) # make sure that tspan_fault has end time >= end time here
include("plotexample_jr.jl")

# simulating a frequency perturbation at node 1
fault1 = ChangeInitialConditions(node = "bus1", var = :ω, f = Inc(0.0));
solution1 = simulate(fault1, powergrid, operationpoint, timespan);
plot1 = create_plot(solution1)

fault1 = ChangeInitialConditions(node = "bus1", var = :ω, f = Inc(0.0));
solution1 = simulate(fault1, powergrid, imported_operationpoint, timespan);
plot1 = create_plot(solution1)