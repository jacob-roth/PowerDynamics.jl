using PowerDynamics: read_powergrid, Json, Inc, find_operationpoint, ChangeInitialConditions, LineFailure, PowerPerturbation, simulate

##
## 118bus with swingeq
##

powergrid = read_powergrid(joinpath(@__DIR__, "grid118.json"), Json)
# operationpoint = find_operationpoint(powergrid, sol_method = :rootfind) # solve_powerflow=true,
# operationpoint = find_operationpoint(powergrid, solve_powerflow = true, sol_method = :rootfind) # solve_powerflow=true,  
# operationpoint = find_operationpoint(powergrid, sol_method = :nlsolve) # solve_powerflow=true, 
operationpoint = find_operationpoint(powergrid, sol_method = :dynamic) # solve_powerflow=true, 

op = import_own_operatingpoint(powergrid, "/Users/jakeroth/git/PowerDynamics.jl/examples/pglib_118bus/")
OP = State(powergrid, op) # convert to `operationpoint` type object

timespan = (0.0, 25) # make sure that tspan_fault has end time >= end time here
include("plotexample_jr.jl")

# simulating a frequency perturbation at node 1
fault1 = ChangeInitialConditions(node = "bus1", var = :ω, f = Inc(0.0));
# solution1 = simulate(fault1, powergrid, OP, timespan); # not working...
solution1 = simulate(fault1, powergrid, operationpoint, timespan);
plot1 = create_plot(solution1)

# simulating a tripped line at line 2
fault2 = LineFailure(line_name = "branch2", tspan_fault = (1.0, 100.0));
solution2 = simulate(fault2, powergrid, operationpoint, timespan);
plot2 = create_plot(solution2)

# confirming that LineFailures object works
fault3 = LineFailures(line_names = ["branch2"], tspan_fault = (1.0, 100.0));
solution3 = simulate(fault3, powergrid, operationpoint, timespan);
plot3 = create_plot(solution3)

# simulating two tripped lines
fault4 = LineFailures(line_names = ["branch2", "branch5"], tspan_fault = (1.0, 100.0));
solution4 = simulate(fault4, powergrid, operationpoint, timespan);
plot4 = create_plot(solution4)


# # simulating a load drop at node 5
# fault3 = PowerPerturbation(node = "bus5", fault_power = 0.0, tspan_fault = (1.0, 5.0), var = :P);
# solution3 = simulate(fault3, powergrid, operationpoint, timespan);
# plot3 = create_plot(solution3)