using PowerDynamics: read_powergrid, Json, Inc, find_operationpoint, ChangeInitialConditions, LineFailure, PowerPerturbation, simulate

##
## 118bus with swingeq
##

powergrid = read_powergrid(joinpath(@__DIR__, "grid118.json"), Json)
operationpoint = find_operationpoint(powergrid, sol_method = :rootfind) # solve_powerflow=true,
operationpoint = find_operationpoint(powergrid, solve_powerflow = true, sol_method = :rootfind) # solve_powerflow=true,  
operationpoint = find_operationpoint(powergrid, sol_method = :nlsolve) # solve_powerflow=true, 
operationpoint = find_operationpoint(powergrid, sol_method = :dynamic) # solve_powerflow=true, 

timespan = (0.0, 5.0)
include("plotexample_jr.jl")
# include("plotexample.jl")

# simulating a frequency perturbation at node 1
fault1 = ChangeInitialConditions(node = "bus1", var = :ω, f = Inc(0.2));
solution1 = simulate(fault1, powergrid, operationpoint, timespan);
plot1 = create_plot(solution1)

# simulating a tripped line between node 1 and 5
fault2 = LineFailure(line_name = "branch2", tspan_fault = (1.0, 5.0));
solution2 = simulate(fault2, powergrid, operationpoint, timespan);
plot2 = create_plot(solution2)

# simulating a load drop at node 5
fault3 = PowerPerturbation(node = "bus5", fault_power = 0.0, tspan_fault = (1.0, 5.0), var = :P);
solution3 = simulate(fault3, powergrid, operationpoint, timespan);
plot3 = create_plot(solution3)

##
## 118bus with swingeq
##
powergrid = read_powergrid(joinpath(@__DIR__, "ieee118_v3.json"), Json)
operationpoint1 = find_operationpoint(powergrid, solve_powerflow = true, sol_method = :rootfind)
operationpoint2 = find_operationpoint(powergrid, solve_powerflow = true, sol_method = :nlsolve)
operationpoint3 = find_operationpoint(powergrid, solve_powerflow = true, sol_method = :dynamic)

timespan = (0.0, 50.0)
include("plotexample_jr.jl")

# simulating a frequency perturbation at node 1
fault1 = ChangeInitialConditions(node = 1, var = :ω, f = Inc(0.0));
solution1 = simulate(fault1, powergrid, operationpoint, timespan);
plot1 = create_plot(solution1)

# simulating a tripped line between node 1 and 5
fault2 = LineFailure(line_name = "branch2", tspan_fault = (1.0, 5.0));
solution2 = simulate(fault2, powergrid, operationpoint, timespan);
plot2 = create_plot(solution2)

# simulating a load drop at node 5
fault3 = PowerPerturbation(node = "bus5", fault_power = 0.0, tspan_fault = (1.0, 5.0), var = :P);
solution3 = simulate(fault3, powergrid, operationpoint, timespan);
plot3 = create_plot(solution3)
