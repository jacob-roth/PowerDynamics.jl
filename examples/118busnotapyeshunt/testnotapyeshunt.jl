include("../../src/PowerDynamics.jl") # comment out `module` part
using MPCCases
include("/Users/jakeroth/git/OPF/src/OPF.jl") # comment out `module` part; this will ERROR, but ok
include("/Users/jakeroth/git/OPF/src/OPF2PD.jl") # comment out `module` part
include("/Users/jakeroth/git/OPF/src/util.jl") # comment out `module` part
include("/Users/jakeroth/git/OPF/src/util_coupled.jl") # comment out `module` part
include("/Users/jakeroth/git/OPF/src/powergrid_utils.jl") # comment out `module` part
include("/Users/jakeroth/git/OPF/src/pfe.jl") # comment out `module` part
using JSON
using PowerDynamics: read_powergrid, Json, Inc, find_operationpoint, ChangeInitialConditions, LineFailure, PowerPerturbation, simulate

##
## 118bus with swingeq
##

powergrid = read_powergrid(joinpath(@__DIR__, "118notapyeshunt.json"), Json)
operationpoint = find_operationpoint(powergrid, sol_method = :dynamic) # solve_powerflow=true, 

operatingdata_path = "/Users/jakeroth/git/PowerDynamics.jl/examples/118busnotapyeshunt/"
imported_vec = import_own_operatingpoint(powergrid, operatingdata_path)
imported_operationpoint = State(powergrid, imported_vec)
norm(imported_operationpoint.vec - operationpoint.vec)

norm(err_from_root(imported_operationpoint)) # 2.942147413731256
norm(err_from_root(operationpoint)) # 2.2147484453872968e-8

timespan = (0.0, 25) # make sure that tspan_fault has end time >= end time here
include("plotexample_jr.jl")

# simulating a frequency perturbation at node 1
fault1 = ChangeInitialConditions(node = "bus1", var = :ω, f = Inc(0.0));
solution1 = simulate(fault1, powergrid, operationpoint, timespan);
plot1 = create_plot(solution1)

fault1 = ChangeInitialConditions(node = "bus1", var = :ω, f = Inc(0.0));
solution1 = simulate(fault1, powergrid, imported_operationpoint, timespan);
plot1 = create_plot(solution1)