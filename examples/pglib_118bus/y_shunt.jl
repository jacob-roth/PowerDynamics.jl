opfmodeldata = get_opfmodeldata(opfdata, case_options)
Y = opfmodeldata[:Y]

casefile_path = "/"
branch_file = readdlm(casefile_path * "mpc_lowdamp_pgliblimits.branch")
branches = Int.(branch_file[:,1:2])

num_branches = size(branches,1)

y_shunt_arr = zeros(Complex, num_branches, num_branches)

for branch_idx in 1:num_branches
    from_bus = branches[branch_idx,1]
    to_bus = branches[branch_idx,2]
    lines = opfmodeldata[:lines]
    buses = opfmodeldata[:buses]
    baseMVA = opfmodeldata[:baseMVA]

    line_b = sum(lines[(lines.from .== from_bus) .& (lines.to .== to_bus)].b)
    from_Bs = first(buses[buses.bus_i .== from_bus].Bs)
    to_Bs = first(buses[buses.bus_i .== to_bus].Bs)

    y_shunt_km = im*(line_b / 2 + (from_Bs / baseMVA))
    y_shunt_mk = im*(line_b / 2 + (to_Bs / baseMVA))
    
    y_shunt_arr[from_bus,to_bus] = y_shunt_km
    y_shunt_arr[to_bus,from_bus] = y_shunt_mk
end

open("/mpc_lowdamp_pgliblimits.y_shunt_arr", "w") do io
    writedlm(io, y_shunt_arr)
end

