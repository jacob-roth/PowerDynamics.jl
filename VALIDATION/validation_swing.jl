import Pkg; Pkg.activate("."); Pkg.instantiate()
using JuMP, Ipopt
using MPCCases
using JSON

include("../src/PowerDynamics.jl") # comment out `module` part
include("/Users/jakeroth/git/OPF/src/OPF.jl") # comment out `module` part; this will ERROR on @constraintref, but ok
include("/Users/jakeroth/git/OPF/src/OPF2PD.jl") # comment out `module` part
include("/Users/jakeroth/git/OPF/src/util.jl") # comment out `module` part
include("/Users/jakeroth/git/OPF/src/util_coupled.jl") # comment out `module` part
include("/Users/jakeroth/git/OPF/src/powergrid_utils.jl") # comment out `module` part
include("/Users/jakeroth/git/OPF/src/pfe.jl") # comment out `module` part
include("/Users/jakeroth/git/OPF/src/powerfloweq.jl")
include("/Users/jakeroth/git/OPF/src/default.jl")
include("/Users/jakeroth/git/OPF/src/util.jl") # comment out `module` part
using PowerDynamics: read_powergrid, Json, Inc, find_operationpoint, ChangeInitialConditions, LineFailure, PowerPerturbation, simulate

include("opt.jl")
include("pd.jl")
include("plot.jl")

##
## get case
##
case_options_base = DefaultOptions();
case_options_base[:lossless]           = true;
case_options_base[:remove_tap]         = true;
case_options_base[:remove_Bshunt]      = false;
case_options_base[:solveOPF]           = false;
case_options_base[:constr_limit_scale] = 1.05;
case_options_base[:current_rating]     = true;
case_options_base[:op_model]           = :acopf;
case_options_base[:shed_load]          = false;
case_path = "/Users/jakeroth/Desktop/planning-large-deviation/data/cases/118-files/"
case_name = "mpc_lowdamp_pgliblimits"
casedata = load_case(case_name, case_path; other=true);
opfdata  = casedata.opf;
case_options = deepcopy(case_options_base);
opfdata   = casedata.opf;
physdata  = casedata.phys;
phys = physDefault;
case_options[:remove_Bshunt] = false
case_options[:nonneg_Bshunt] = false
case_options[:nobus_Bshunt] = true

opfmodeldata = get_opfmodeldata(opfdata, case_options);
opfmodeldata[:Y] = imag.(opfmodeldata[:Y])
nline = length(opfmodeldata[:lines])

## albert's operating data
operating_data_path = "/Users/jakeroth/git/PowerDynamics.jl/src/simple_kmc_validation/1_00__0__0__acopf__1e_09/"
optimal_values = Dict()
optimal_values[:Pnet] = readdlm(operating_data_path * "Pnet.csv")
optimal_values[:Qnet] = readdlm(operating_data_path * "Qnet.csv")
optimal_values[:Vm] = readdlm(operating_data_path * "Vm.csv")
optimal_values[:Va] = readdlm(operating_data_path * "Va.csv")

##
## run simple failures
##
op_kmc_dict, status = simple_energy_model(opfmodeldata, optimal_values)
oklinenames, oklinemask = get_ok_lines(powergrid0)
oklines = parse.(Int,[replace(k,"branch" => "" ) for k in oklinenames])
op_kmc_dict_ls = []
status_ls = []
for l in oklines
    op_l, status_l = simple_energy_model_failure(opfdata, case_options, optimal_values, l)
    push!(op_kmc_dict_ls, op_l)
    push!(status_ls, status_l)
end

id2idx = Dict()
for l in oklines
    id2idx[l] = findfirst(oklines .== l)
end
@assert(id2idx[13] == 3)

##
## convert to PD and get failure times
##
grid_file = "/Users/jakeroth/git/PowerDynamics.jl/VALIDATION/ieee118_fpacopf_1e09_swing.json"
powergrid0 = read_powergrid(grid_file, Json)
op_kmc_vec = import_own_operatingpoint(powergrid0, op_kmc_dict)
op_kmc = State(powergrid0, op_kmc_vec)
@assert(norm(err_from_root(op_kmc)) < 5e-5)

op_pd = find_operationpoint(powergrid0, sol_method=:dynamic)
@assert(norm(err_from_root(op_pd)) < 5e-5)

global op_pds = Array{Any}(undef, length(oklines),1)
global tstable_op_pds = Array{Any}(undef, length(oklines), 1)
tol = 1e-5
for l in oklines
    println("line=$l, idx=$(id2idx[l])")
    fault_l = LineFailure(line_name="branch$(l)", tspan_fault=(1,Inf))
    try
        global sol_pd = simulate(fault_l, powergrid0, op_pd, (0,1e3));
        global tstable_op_pd = find_stable_time(sol_pd, powergrid0, fault_l, tol)
        global op_pds[id2idx[l]] = sol_pd.dqsol[end];
        global tstable_op_pds[id2idx[l]] = tstable_op_pd;
    catch e
        println(e)
        global sol_pd = [-Inf]
        global tstable_op_pd = -Inf
    finally
        global op_pds[id2idx[l]] = sol_pd
        global tstable_op_pds[id2idx[l]] = tstable_op_pd
    end
end

writedlm("data/op_pds_swing.csv", op_pds)
writedlm("data/tstable_op_pds_swing.csv", tstable_op_pds)

# op_pds = readdlm("data/op_pds_swing.csv")
# tstable_op_pds = readdlm("data/tstable_op_pds_swing.csv")

pg0 = deepcopy(powergrid0)
unstable = findall(vec(tstable_op_pds) .< 0) # SAME AS KMC RIGHT NOW...
stable = findall(vec(tstable_op_pds) .> 0) # SAME AS KMC RIGHT NOW...
op_stability_kmc = []
# op_stability_pd = []
for l in oklines
    idx = id2idx[l]
    println("line $l")
    fault_l = LineFailure(line_name="branch$(l)", tspan_fault=(1,Inf))
    pg = deepcopy(pg0)
    grid_l = fault_l(pg)
    # op_pd_l_vec = op_pds[idx,:]
    op_kmc_l_vec = import_own_operatingpoint(grid_l, op_kmc_dict_ls[idx])
    # op_pd_l_state = State(grid_l, op_pd_l_vec)
    op_kmc_l_state = State(grid_l, op_kmc_l_vec)
    err_kmc_vec = err_from_root(op_kmc_l_state)
    # err_pd_vec = err_from_root(op_pd_l_state)
    err_kmc = norm(err_kmc_vec)
    # err_pd = norm(err_pd_vec)
    push!(op_stability_kmc, err_kmc)
    # push!(op_stability_pd, err_pd)
end

writedlm("data/op_stability_kmc_swing.csv", op_stability_kmc)


##
## stats
##
maximum(tstable_op_pds)
mean(tstable_op_pds[tstable_op_pds.>0])

powergrid = powergrid0
PD_line_names = collect(keys(powergrid.lines))
PD_line_IDs = parse.(Int,[replace(k,"branch" => "" ) for k in PD_line_names])

err_IDs = findall(vec(tstable_op_pds).<0)
nonerr_IDs = findall(vec(tstable_op_pds).>0)
stable_line_IDs = nonerr_IDs
unstable_line_IDs = err_IDs

## pd totals
n_ok_PD_lines = length(oklines)
n_stable_lines = length(stable_line_IDs)

## kmc stable/unstable for PD stable
n_KMC_unstable_PD_stable = sum(abs.(op_stability_kmc[stable_line_IDs]).>1e-5)
n_KMC_stable_PD_stable = sum(abs.(op_stability_kmc[stable_line_IDs]).<=1e-5)

## kmc stable/unstable for PD unstable
n_KMC_stable_PD_unstable = sum(abs.(op_stability_kmc[unstable_line_IDs]).>1e-5)
n_KMC_unstable_PD_unstable = sum(abs.(op_stability_kmc[unstable_line_IDs]).<=1e-5)


op_pds = Array{Any}(undef, 5)
tstable_op_pds = Array{Any}(undef, 5)
for l in 1:2;
    if l==2; throw("er"); end;
    try
        global sol_pd2 = 1
    catch e
        println(e)
        global sol_pd2 = 2
    finally
        op_pds[l] = sol_pd2
    end
end