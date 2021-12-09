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
op_kmc_dict_ls = []
status_ls = []
for l in 1:nline
    op_l, status_l = simple_energy_model_failure(opfdata, case_options, optimal_values, l)
    push!(op_kmc_dict_ls, op_l)
    push!(status_ls, status_l)
end

##
## convert to PD and get failure times
##
grid_file = "/Users/jakeroth/git/PowerDynamics.jl/examples/pglib_118bus/ieee118_fpacopf_1e09.json"
powergrid0 = read_powergrid(grid_file, Json)
op_kmc_vec = import_own_operatingpoint(powergrid0, op_kmc_dict)
op_kmc = State(powergrid0, op_kmc_vec)
@assert(norm(err_from_root(op_kmc)) < 5e-5)

op_pd = find_operationpoint(powergrid0, sol_method=:dynamic)
@assert(norm(err_from_root(op_pd)) < 5e-5)

op_pds = []
op_kmcs = []
tstable_op_pds = []
tstable_op_kmcs = []
op_diffs = []
tol = 5e-5
for l in 1:nline
    println(l)
    fault_l = LineFailure(line_name="branch$(l)", tspan_fault=(1,Inf))
    try
        sol_kmc = simulate(fault_l, powergrid0, op_kmc, (0,1e9));
        sol_pd = simulate(fault_l, powergrid0, op_pd, (0,1e9));
    catch e
        println(e)
    finally
        # create_plot(sol_kmc, powergrid0)
        # create_plot(sol_pd, powergrid0)
        tstable_op_pd = find_stable_time(sol_pd, powergrid0, fault_l, tol)
        tstable_op_kmc = find_stable_time(sol_kmc, powergrid0, fault_l, tol)
        op_diff = norm(sol_pd.dqsol[end] .- sol_kmc.dqsol[end])
        push!(op_pds, sol_pd.dqsol[end])
        push!(op_kmcs, sol_kmc.dqsol[end])
        push!(tstable_op_pds, tstable_op_pd)
        push!(tstable_op_kmcs, tstable_op_kmc)
        push!(op_diffs, op_diff)
    end
end

writedlm("data/op_pds.csv", op_pds)
writedlm("data/op_kmcs.csv", op_kmcs)
writedlm("data/tstable_op_pds.csv", tstable_op_pds)
writedlm("data/tstable_op_kmcs.csv", tstable_op_kmcs)
writedlm("data/op_diffs.csv", op_diffs)

op_pds = readdlm("data/op_pds.csv")
op_kmcs = readdlm("data/op_kmcs.csv")
tstable_op_pds = readdlm("data/tstable_op_pds.csv")
tstable_op_kmcs = readdlm("data/tstable_op_kmcs.csv")
op_diffs = readdlm("data/op_diffs.csv")

function get_ok_lines(powergrid)
    ok = []
    for lkey in keys(powergrid.lines)
        f = powergrid.lines[lkey].from
        t = powergrid.lines[lkey].to
        if typeof(powergrid.nodes[f]) == SwingEqLVS || typeof(powergrid.nodes[t]) == SwingEqLVS
            continue
        end
        push!(ok, lkey)
    end
    return ok
end

fig, ax = subplots(figsize=(4,4))
nok = sum(tstable_op_pds .> 0)
ax.cla()
ax.hist(tstable_op_pds[tstable_op_pds .> 0],   alpha=0.5, label="PD")
ax.hist(tstable_op_kmcs[tstable_op_kmcs .> 0], alpha=0.5,  label="KMC")
ax.set_title("Time to Stability: $(nok) / $(nline) lines")
ax.legend()
ax.set_xlabel("time (s)")
ax.set_ylabel("number of lines")
fig.savefig("time_to_stability.pdf")

ax.cla()
ax.hist(op_diffs)

op_kmc_dict_ls