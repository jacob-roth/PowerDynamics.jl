import Pkg;
Pkg.activate(".");
Pkg.instantiate();
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
# using PowerDynamics: read_powergrid, Json, Inc, find_operationpoint, ChangeInitialConditions, LineFailure, PowerPerturbation, simulate

include("opt.jl")
include("pd.jl")
include("plot.jl")

##
## get case
##
case_options_base = DefaultOptions();
case_options_base[:lossless] = true;
case_options_base[:remove_tap] = true;
case_options_base[:remove_Bshunt] = false;
case_options_base[:solveOPF] = false;
case_options_base[:constr_limit_scale] = 1.05;
case_options_base[:current_rating] = true;
case_options_base[:op_model] = :acopf;
case_options_base[:shed_load] = false;
case_path = "/Users/jakeroth/Desktop/planning-large-deviation/data/cases/118-files/"
case_name = "mpc_lowdamp_pgliblimits"
casedata = load_case(case_name, case_path; other = true);
opfdata = casedata.opf;
case_options = deepcopy(case_options_base);
opfdata = casedata.opf;
physdata = casedata.phys;
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
grid_file = "/Users/jakeroth/git/PowerDynamics.jl/VALIDATION/ieee118_fpacopf_1e09_swing.json"
powergrid0 = read_powergrid(grid_file, Json)

op_kmc_dict, status = simple_energy_model(opfmodeldata, optimal_values)
oklinenames, oklinemask = get_ok_lines(powergrid0)
oklines = parse.(Int, [replace(k, "branch" => "") for k in oklinenames])
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
op_kmc_vec = import_own_operatingpoint(powergrid0, op_kmc_dict)
op_kmc = State(powergrid0, op_kmc_vec)
@assert(norm(err_from_root(op_kmc)) < 5e-5)

op_pd = find_operationpoint(powergrid0, sol_method = :dynamic)
@assert(norm(err_from_root(op_pd)) < 5e-5)

global op_kmc2pds = Array{Any}(undef, length(oklines))
global tstable_op_kmc2pds = Array{Any}(undef, length(oklines))
tol = 1e-5
for l in oklines
    println("line=$l, idx=$(id2idx[l])")
    fault_l = LineFailure(line_name = "branch$(l)", tspan_fault = (1, Inf))
    try
        global sol_kmc2pd = simulate(fault_l, powergrid0, op_kmc, (0, 1e3))
        global tstable_op_kmc2pd = find_stable_time(sol_kmc2pd, powergrid0, fault_l, tol)
        global op_kmc2pds[id2idx[l]] = sol_kmc2pd.dqsol[end]
        global tstable_op_kmc2pds[id2idx[l]] = tstable_op_kmc2pd
    catch e
        println(e)
        global sol_pd = -Inf
        global tstable_op_pd = -Inf
    finally
        global op_kmc2pds[id2idx[l]] = sol_kmc2pd.dqsol[end]
        global tstable_op_kmc2pds[id2idx[l]] = tstable_op_kmc2pd
    end
end

writedlm("data/op_kmc2pds_swing.csv", op_kmc2pds)
writedlm("data/tstable_op_kmc2pds_swing.csv", tstable_op_kmc2pds)

##
## compare operating points
##

op_kmc2pds = readdlm("data/op_kmc2pds_swing.csv")
tstable_op_kmc2pds = readdlm("data/tstable_op_kmc2pds_swing.csv")

pg0 = deepcopy(powergrid0)
unstable = findall(vec(tstable_op_kmc2pds) .< 0)
stable = findall(vec(tstable_op_kmc2pds) .> 0)
op_stability_kmc2pd = []
op_stability_kmc2kmc = []
op_diffs = []
for l in oklines
    idx = id2idx[l]
    println("line $l")
    fault_l = LineFailure(line_name = "branch$(l)", tspan_fault = (1, Inf))
    pg = deepcopy(pg0)
    grid_l = fault_l(pg)
    kmc2pd_v = op_kmc2pds[idx, :] # solve PD from line failures for op point
    kmc2kmc_v = import_own_operatingpoint(grid_l, op_kmc_dict_ls[idx]) # solve EM (energy model) for line failures (as per kmc) for op point
    kmc2kmc_state = State(grid_l, kmc2kmc_v)
    kmc2pd_state = State(grid_l, kmc2pd_v)

    err_kmc2kmc_vec = err_from_root(kmc2kmc_state)
    err_kmc2pd_vec = err_from_root(kmc2pd_state)
    err_kmc2kmc = norm(err_kmc2kmc_vec)
    err_kmc2pd = norm(err_kmc2pd_vec)
    push!(op_stability_kmc2kmc, err_kmc2kmc) # is the point solved by EM PD-stable?
    push!(op_stability_kmc2pd, err_kmc2pd) # is the point solved by PD PD-stable?
    push!(op_diffs, norm(kmc2kmc_v - kmc2pd_v)) # are the points the same?
end

writedlm("data/op_stability_kmc2pd_swing.csv", op_stability_kmc2pd)
writedlm("data/op_stability_kmc2kmc_swing.csv", op_stability_kmc2kmc)

##
## stats: IV.D.3 paragraph 2
##

# energy model
stable_lines_em = oklines[abs.(vec(op_stability_kmc2kmc)).<1e-5] # 127 lines
unstable_lines_em = oklines[abs.(vec(op_stability_kmc2kmc)).>1e-5] # 7 lines

# power dynamics
stable_lines_pd = oklines[vec(tstable_op_kmc2pds).>0] # 129 lines
unstable_lines_pd = oklines[vec(tstable_op_kmc2pds).<0] # 5 lines

# both
stable_lines_both = oklines[(abs.(vec(op_stability_kmc2kmc)).<1e-5).*(vec(tstable_op_kmc2pds).>0)] # 126 lines
unstable_lines_both = oklines[(abs.(vec(op_stability_kmc2kmc)).>1e-5).*(vec(tstable_op_kmc2pds).<0)] # 4 lines
stable_em_unstable_pd = oklines[(abs.(vec(op_stability_kmc2kmc)).<1e-5).*(vec(tstable_op_kmc2pds).<0)] # 1 lines
unstable_em_stable_pd = oklines[(abs.(vec(op_stability_kmc2kmc)).>1e-5).*(vec(tstable_op_kmc2pds).>0)] # 3 lines

# true/false postive/negative rates
em_tpr_stable = length(stable_lines_both) / length(stable_lines_pd) # 98% = 126 / 129
em_fnr_stable = length(unstable_em_stable_pd) / length(stable_lines_pd) # 2% = 3 / 129
em_fpr_unstable = length(stable_em_unstable_pd) / length(unstable_lines_pd) # 20% = 1 / 5
em_tpr_unstable = length(unstable_lines_both) / length(unstable_lines_pd) # 80% = 4 / 5

# overall stability pct
stable_overall_pct_pd = length(stable_lines_pd) / length(oklines) # 96% = 129 / 134; PD-stable
stable_overall_pct_em = length(stable_lines_em) / length(oklines) # 95% = 127 / 134; EM-stable

# timing stats (for PD-stable)
maximum(tstable_op_kmc2pds) # 34.6s
mean(tstable_op_kmc2pds[tstable_op_kmc2pds.>0]) # 27.4s
std(tstable_op_kmc2pds[tstable_op_kmc2pds.>0]) # 3.1s

# when didn't they find the same operating point?
mask = abs.(op_diffs) .> 1e-5
diff_lines__stable_em_stable_pd = oklines[mask.*(abs.(vec(op_stability_kmc2kmc)).<1e-5).*(vec(tstable_op_kmc2pds).>0)]
diff_lines__unstable_em_stable_pd = oklines[mask.*(abs.(vec(op_stability_kmc2kmc)).>1e-5).*(vec(tstable_op_kmc2pds).>0)]
diff_lines__stable_em_unstable_pd = oklines[mask.*(abs.(vec(op_stability_kmc2kmc)).<1e-5).*(vec(tstable_op_kmc2pds).<0)]
diff_lines__unstable_em_unstable_pd = oklines[mask.*(abs.(vec(op_stability_kmc2kmc)).>1e-5).*(vec(tstable_op_kmc2pds).<0)]
