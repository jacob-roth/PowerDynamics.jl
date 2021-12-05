##
## 30bus case: activate / run from the `OPF/tests` environment
##
using JuMP, Ipopt
using MPCCases
using JSON
include("/Users/jakeroth/git/OPF/src/OPF.jl") # uncommented module
path = "/Users/jakeroth/Desktop/planning-large-deviation/data/cases/30-files/"
tol = 1e-6
plotting = false
temp = 1e-4
damping = 1.0
case_name = "mpc-lowdamp"

options = DefaultOptions()
options[:current_rating] = true
options[:lossless]       = true
options[:remove_Bshunt]  = true
options[:remove_tap]     = true
options[:shed_load]      = false
options[:print_level]    = 1
options[:constr_limit_scale] = 1.05
options[:pw_angle_limits]= false
options[:slack0]         = true

casedata = load_case(case_name, path; other=true)
opfdata  = casedata.opf
opfmodel = acopf_model(opfdata, options)
opfmodel_acopf = acopf_solve(opfmodel, opfdata)
opfmodeldata = get_opfmodeldata(opfdata, options)
opfmodeldata[:Y] = imag.(opfmodeldata[:Y])
solution = get_optimal_values(opfmodel_acopf.m, opfmodeldata)
include("/Users/jakeroth/git/OPF/src/OPF2PD.jl")
line_type = "PiModelLine"

opf2pd("/Users/jakeroth/git/PowerDynamics.jl/examples/30bus/30bus.json", solution, opfmodeldata, line_type, physDefault)
writedlm("/Users/jakeroth/git/PowerDynamics.jl/examples/30bus/Pnet.csv", solution[:Pnet])
writedlm("/Users/jakeroth/git/PowerDynamics.jl/examples/30bus/Qnet.csv", solution[:Qnet])
writedlm("/Users/jakeroth/git/PowerDynamics.jl/examples/30bus/Vm.csv", solution[:Vm])
writedlm("/Users/jakeroth/git/PowerDynamics.jl/examples/30bus/Va.csv", solution[:Va])
writedlm("/Users/jakeroth/git/PowerDynamics.jl/examples/30bus/Y.csv", opfmodeldata[:Y])