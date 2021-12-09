##
## environment
##
import Pkg; Pkg.activate("."); Pkg.instantiate()
include("../../src/PowerDynamics.jl") # comment out `module` part
using MPCCases
using JSON
include("/Users/jakeroth/git/OPF/src/OPF.jl") # comment out `module` part; this will ERROR on @constraintref, but ok
include("/Users/jakeroth/git/OPF/src/OPF2PD.jl") # comment out `module` part
include("/Users/jakeroth/git/OPF/src/util.jl") # comment out `module` part
include("/Users/jakeroth/git/OPF/src/util_coupled.jl") # comment out `module` part
include("/Users/jakeroth/git/OPF/src/powergrid_utils.jl") # comment out `module` part
include("/Users/jakeroth/git/OPF/src/pfe.jl") # comment out `module` part
using PowerDynamics: read_powergrid, Json, Inc, find_operationpoint, ChangeInitialConditions, LineFailure, PowerPerturbation, simulate
##
## environment
##

##
## failure info
##
failure_idx = 1
seq_idx = 1
num_seqs = 10
gridname = "seq_$(seq_idx)_failure_$(failure_idx).json"
##
## failure info
##


##
## case data
##
# operatingdata_path = "/Users/Albert/Personal/Argonne/cascades/kmcdata/118/paper/review_2/seq_10/v3/"
operating_data_path = "/Users/jakeroth/gdrive/review2/seq_10_v3/"
kmc_dict = load_kmc_dict(operating_data_path, failure_idx)

# fileout = "/Users/Albert/Personal/Argonne/cascades/casedata/118/paper_experiments/seq_$(seq_idx)_failure_$(failure_idx).json"
fileout = "/Users/jakeroth/git/PowerDynamics.jl/src/kmc_validation/kmc_grid_data/$(gridname)"
line_type = "PiModelLine"
bus_type = "SwingEqLVS"
nobus_Bshunt = true

casedata_path = "/Users/jakeroth/git/PowerDynamics.jl/examples/118bus_nobus_Bshunt/"
case_options_base = DefaultOptions();
case_options_base[:lossless]           = true;
case_options_base[:remove_tap]         = true;
case_options_base[:remove_Bshunt]      = false;
case_options_base[:solveOPF]           = false;
case_options_base[:constr_limit_scale] = 1.05;
case_options_base[:current_rating]     = true;
case_options_base[:op_model]           = :acopf;
case_options_base[:shed_load]          = false;

casepath = "/Users/jakeroth/Desktop/planning-large-deviation/data/cases/118-files/"
casename = "mpc_lowdamp_pgliblimits"
casedata = load_case(casename, casepath; other=true);
opfdata  = casedata.opf;
case_options = deepcopy(case_options_base);
opfdata   = casedata.opf;
physdata  = casedata.phys;

case_options[:remove_Bshunt] = false
case_options[:nonneg_Bshunt] = false
case_options[:nobus_Bshunt] = true

opfmd = get_opfmodeldata(opfdata, case_options);
kmc2pd(fileout, opfmd, kmc_dict, line_type, bus_type, seq_idx, nobus_Bshunt)
##
## case data
##


##
## validation experiment
##
degraded_powergrid = read_powergrid(fileout, Json)
op_kmc_vec = import_own_operatingpoint(degraded_powergrid, kmc_dict)
op_kmc = State(degraded_powergrid, op_kmc_vec)
err_from_root(op_kmc)

fault1 = ChangeInitialConditions(node = "bus1", var = :ω, f = Inc(0.0));
solution1 = simulate(fault1, degraded_powergrid, op_kmc, (0,10.0));
plot1 = create_plot(solution1)

##
## validation experiment
##