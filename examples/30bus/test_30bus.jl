include("../../src/PowerDynamics.jl") # comment out `module` part
include("/Users/jakeroth/git/OPF/src/OPF2PD.jl") # comment out `module` part
include("/Users/jakeroth/git/OPF/src/util.jl") # comment out `module` part
include("/Users/jakeroth/git/OPF/src/util_coupled.jl") # comment out `module` part
include("/Users/jakeroth/git/OPF/src/powergrid_utils.jl") # comment out `module` part
include("/Users/jakeroth/git/OPF/src/pfe.jl") # comment out `module` part
using JSON

##
## load power grid and get operating point
##

operatingdata_path = "/Users/jakeroth/git/PowerDynamics.jl/examples/30bus/"
powergrid = read_powergrid(joinpath(@__DIR__, "30bus.json"), Json)
operationpoint = find_operationpoint(powergrid, sol_method = :dynamic)

imported_vec = import_own_operatingpoint(powergrid, operatingdata_path)
imported_operationpoint = State(powergrid, imported_vec)
norm(imported_operationpoint.vec - operationpoint.vec)

omega_import, ur_import, ui_import = get_components(imported_operationpoint)
omega, ur, ui = get_components(operationpoint)
v = ur .+ im.*ui
v_import = ur_import .+ im.*ui_import
vm, va = abs.(v), angle.(v)
vm_import, va_import = abs.(v_import), angle.(v_import)

norm(vm-vm_import)
norm(va-va_import)

Pnet = readdlm(operatingdata_path * "Pnet.csv")
Qnet = readdlm(operatingdata_path * "Qnet.csv")
P, Q = powerfloweqs(vm,va,opfmodeldata)
P_import, Q_import = powerfloweqs(vm_import, va_import, opfmodeldata)

P - P_import .> 1e-5
Q - Q_import .> 1e-5

##
## simulate
##

timespan = (0,100.0)
fault1 = ChangeInitialConditions(node = "bus2", var = :ω, f = Inc(0.0));
solution1 = simulate(fault1, powergrid, operationpoint, timespan);
plot1 = create_plot(solution1)

timespan = (0,100.0)
fault1 = ChangeInitialConditions(node = "bus2", var = :ω, f = Inc(0.0));
solution1 = simulate(fault1, powergrid, imported_operationpoint, timespan);
plot1 = create_plot(solution1)
