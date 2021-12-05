include("../../src/PowerDynamics.jl") # comment out `module` part
include("/Users/jakeroth/git/OPF/src/OPF2PD.jl") # comment out `module` part
include("/Users/jakeroth/git/OPF/src/util.jl") # comment out `module` part
include("/Users/jakeroth/git/OPF/src/util_coupled.jl") # comment out `module` part
include("/Users/jakeroth/git/OPF/src/powergrid_utils.jl") # comment out `module` part
include("/Users/jakeroth/git/OPF/src/pfe.jl") # comment out `module` part
using JSON

##
## get power grid data
##
operatingdata_path = "/Users/jakeroth/git/PowerDynamics.jl/examples/pglib_118bus/1_00__0__0__acopf__1_05/"
line_type = "PiModelLine"#Tapratio"
line_type = "PiModelLineTapratio"

casepath = "/Users/jakeroth/Desktop/planning-large-deviation/data/cases/118-files/"
casename = "mpc_lowdamp_pgliblimits"
casedata = load_case(casename, casepath; other=true)
opfdata  = casedata.opf

options = DefaultOptions()
options[:emergencylimit] = NaN
options[:ratelimit]      = NaN
options[:current_rating] = true
options[:lossless]       = true
options[:remove_Bshunt]  = false
options[:remove_tap]     = false
options[:shed_load]      = false
options[:print_level]    = 1
options[:constr_limit_scale] = 1.05
options[:pw_angle_limits]= false
options[:slack0]         = true

opfmodeldata = get_opfmodeldata(opfdata, options)
opfmodeldata[:nbus] = length(opfdata.buses)
opfmodeldata[:Y] = imag.(opfmodeldata[:Y])

# write out
fileout = "/Users/jakeroth/git/PowerDynamics.jl/examples/pglib_118bus/grid118_new.json"
opf2pd(fileout, operatingdata_path, opfmodeldata, line_type)

##
## load power grid and get operating point
##

powergrid = read_powergrid(joinpath(@__DIR__, "grid118_new.json"), Json)
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

lineratios = findall(lines.ratio.!=0)
findall(P - P_import .> 1e-5)