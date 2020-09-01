using Ipopt: Optimizer
using PowerModels: ACPPowerModel, run_pf

function make_generator!(dict::Dict{String,Any}, key_b::Int, node)
    key = length(dict["gen"]) + 1
    (dict["gen"])[string(key)] = Dict{String,Any}()
    ((dict["gen"])[string(key)])["mBase"] = 1
    ((dict["gen"])[string(key)])["gen_bus"] = key_b
    ((dict["gen"])[string(key)])["gen_status"] = 1
    ((dict["gen"])[string(key)])["source_id"] = Any["gen", key]
    ((dict["gen"])[string(key)])["index"] = key

    if isa(node, SlackAlgebraic)
        ((dict["gen"])[string(key)])["vg"] = node.U
        # ((dict["gen"])[string(key)])["pg"] = no
        ((dict["gen"])[string(key)])["pmin"] = 0
        # ((dict["gen"])[string(key)])["pmax"] = 1.1 * node.P
    else
        ((dict["gen"])[string(key)])["pg"] = node.P
        ((dict["gen"])[string(key)])["pmin"] = 0.9 * node.P
        ((dict["gen"])[string(key)])["pmax"] = 1.1 * node.P
        ((dict["gen"])[string(key)])["vg"] = node.V
    end

    ((dict["gen"])[string(key)])["qg"] = 0
    ((dict["gen"])[string(key)])["qmin"] = -1.5
    ((dict["gen"])[string(key)])["qmax"] = 1.5
end

"""
The entries that are present for all busses.
"""
function _make_bus_ac_header(data::Dict{String,Any})
    key_e = length(data["bus"]) + 1
    (data["bus"])[string(key_e)] = Dict{String,Any}()
    bus_dict = (data["bus"])[string(key_e)]
    bus_dict["source_id"] = Any["bus", key_e]
    bus_dict["index"] = key_e
    bus_dict["vmin"] = 0.9
    bus_dict["vmax"] = 1.1
    bus_dict["vm"] = 1
    bus_dict["va"] = 0
    #bus_dict["bus_i"] = key_e #bus_dict["zone"] = 1 #bus_dict["area"] = 1 #bus_dict["base_kv"] = 1
    return bus_dict
end

"""
The default case.
"""
function make_bus_ac!(data::Dict{String,Any}, node)
    bus_dict = _make_bus_ac_header(data)
    bus_dict["bus_type"] = 4
end

function make_bus_ac!(data::Dict{String,Any}, node::PQAlgebraic)
    bus_dict = _make_bus_ac_header(data)
    bus_dict["bus_type"] = 1

    key_l = length(data["load"]) + 1
    (data["load"])[string(key_l)] = Dict{String,Any}()
    dict_l = (data["load"])[string(key_l)]
    dict_l["source_id"] = bus_dict["source_id"]
    dict_l["load_bus"] = bus_dict["index"]
    dict_l["status"] = 1
    dict_l["pd"] = node.P
    dict_l["qd"] = node.Q
    dict_l["index"] = key_l
end

function make_bus_ac!(data::Dict{String,Any}, node::Union{PVAlgebraic,SwingEqLVS})
    bus_dict = _make_bus_ac_header(data)
    bus_dict["bus_type"] = 2
    bus_dict["vm"] = abs(node.V)  # assumed p.u.
    bus_dict["vmin"] = 0.9 * bus_dict["vm"]
    bus_dict["vmax"] = 1.1 * bus_dict["vm"]
    bus_dict["va"] = angle(node.V)
    isa(node, SwingEqLVS) ? make_generator!(data, bus_dict["index"], node) : nothing
end

function make_bus_ac!(data::Dict{String,Any}, node::SlackAlgebraic)
    bus_dict = _make_bus_ac_header(data)
    bus_dict["bus_type"] = 3
    bus_dict["vm"] = abs(node.U)  # assumed p.u.
    bus_dict["vmin"] = 0.9 * bus_dict["vm"]
    bus_dict["vmax"] = 1.1 * bus_dict["vm"]
    bus_dict["va"] = angle(node.U)
    make_generator!(data, bus_dict["index"], node)
end

function _make_branch_ac_header(data::Dict{String,Any})
    key_e = length(data["branch"]) + 1
    (data["branch"])[string(key_e)] = Dict{String,Any}()
    branch_dict = (data["branch"])[string(key_e)]
    branch_dict["source_id"] = Any["branch", key_e]
    branch_dict["index"] = key_e
    branch_dict["rate_a"] = 1
    branch_dict["rate_b"] = 1
    branch_dict["rate_c"] = 1
    branch_dict["c_rating_a"] = 1
    branch_dict["br_status"] = 1
    branch_dict["angmin"] = ang_min
    branch_dict["angmax"] = ang_max

    # default values
    branch_dict["transformer"] = false
    branch_dict["tap"] = 1
    branch_dict["shift"] = 0
    branch_dict["g_fr"] = 0
    branch_dict["b_fr"] = 0
    branch_dict["g_to"] = 0
    branch_dict["b_to"] = 0
    return branch_dict
end

"""
The default, non-specialised method.
"""
function make_branch_ac!(data::Dict{String,Any}, line)
    throw(ArgumentError("Line type $(typeof(line)) does not exist."))
end

function make_branch_ac!(data::Dict{String,Any}, line::StaticLine)
    branch_dict = _make_branch_ac_header(data)
    branch_dict["f_bus"] = line.from
    branch_dict["t_bus"] = line.to
    branch_dict["br_r"] = real(1 / line.Y)
    branch_dict["br_x"] = imag(1 / line.Y)
end

function make_branch_ac!(data::Dict{String,Any}, line::RLLine)
    branch_dict = _make_branch_ac_header(data)
    branch_dict["f_bus"] = line.from
    branch_dict["t_bus"] = line.to
    branch_dict["br_r"] = real(line.R)
    branch_dict["br_x"] = imag(line.L * line.ω0)
end

function make_branch_ac!(data::Dict{String,Any}, line::Transformer)
    branch_dict = _make_branch_ac_header(data)
    branch_dict["f_bus"] = line.from
    branch_dict["t_bus"] = line.to
    branch_dict["transformer"] = true
    branch_dict["tap"] = line.t_ratio
    branch_dict["br_r"] = real(1 / line.Y)
    branch_dict["br_x"] = imag(1 / line.Y)
end

# PiModel is not a line type but a function that the admittance matrix.
# The following will not work:
# function make_branch_ac!(data :: Dict{String, Any}, line::PiModel)
#     _make_branch_ac_header(data)
#     data["transformer"] = true
#     data["tap"] = line.t_km / line.t_mk
#     data["g_fr"] = real(line.y_shunt_km / data["tap"]^2)
#     data["b_fr"] = imag(line.y_shunt_km / data["tap"]^2)
#     data["br_r"] = real(1/line.Y)
#     data["br_x"] = imag(1/line.Y)
#     data["g_to"] = real(line.y_shunt_mk)
#     data["b_to"] = imag(line.y_shunt_mk)
# end

function make_branch_ac!(data::Dict{String,Any}, line::PiModelLine)
    branch_dict = _make_branch_ac_header(data)
    branch_dict["f_bus"] = line.from
    branch_dict["t_bus"] = line.to
    branch_dict["g_fr"] = real(line.y_shunt_km)
    branch_dict["b_fr"] = imag(line.y_shunt_km)
    branch_dict["br_r"] = real(1 / line.Y)
    branch_dict["br_x"] = imag(1 / line.Y)
    branch_dict["g_to"] = real(line.y_shunt_mk)
    branch_dict["b_to"] = imag(line.y_shunt_mk)
end


function power_flow(power_grid::PowerGrid)
    global ang_min, ang_max
    ang_min = deg2rad(60)
    ang_max = deg2rad(-60)

    data = Dict{String,Any}()
    data["source_type"] = "PowerDynamics"
    data["name"] = "network"
    data["source_version"] = "2.3.2"
    data["per_unit"] = true
    data["baseMVA"] = 1
    keywords = [
        "bus"
        "shunt"
        "dcline"
        "storage"
        "switch"
        "load"
        "branch"
        "gen"
    ]

    for keyword in keywords
        data[keyword] = Dict{String,Any}()
    end

    for line in power_grid.lines
        make_branch_ac!(data, line)
    end

    for node in power_grid.nodes
        make_bus_ac!(data, node)
    end

    result = run_pf(data, ACPPowerModel, Optimizer)

    return data, result
end

export power_flow
