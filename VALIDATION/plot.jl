# using Plots
# using PowerDynamics: SwingEq

function create_plot(sol, powergrid)
    generator_indices = findall(bus -> typeof(bus) == SwingEq, powergrid.nodes)
    labels = reshape(generator_indices,(1,length(generator_indices)))

    pl_v = Plots.plot(sol, generator_indices, :v, legend = (0.8, 0.7), ylabel="V [p.u.]",  label=labels)
    pl_p = Plots.plot(sol, generator_indices, :p, legend = (0.8, 0.7), ylabel="p [p.u.]",  label=labels)
    pl_q = Plots.plot(sol, generator_indices, :q, legend = (0.8, 0.7), ylabel="q [p.u.]",  label=labels)
    pl_ω = Plots.plot(sol, generator_indices, :ω, legend = (0.8, 0.7), ylabel="ω [rad/s]", label=labels)

    pl = Plots.plot(pl_ω, pl_v, pl_p, pl_q;
            layout=(2,2),
            size = (1000, 500),
            lw=3,
            xlabel="t[s]")
end

using PyPlot
function create_plot(sol, powergrid)
    generator_indices = findall(bus -> typeof(bus) == SwingEqLVS, powergrid.nodes)
    labels = reshape(generator_indices,(1,length(generator_indices)))
    gids = parse.(Int,[replace(k,"bus" => "" ) for k in generator_indices])
    fig, ax = plt.subplots()

    
    pl_v = ax.plot([x.u[gids][1] for x in sol.dqsol])
    pl_p = ax.plot(sol.dqsol, generator_indices, :p, legend = (0.8, 0.7), ylabel="p [p.u.]",  label=labels)
    pl_q = ax.plot(sol.dqsol, generator_indices, :q, legend = (0.8, 0.7), ylabel="q [p.u.]",  label=labels)
    pl_ω = ax.plot(sol.dqsol, generator_indices, :ω, legend = (0.8, 0.7), ylabel="ω [rad/s]", label=labels)

    pl = Plots.plot(pl_ω, pl_v, pl_p, pl_q;
            layout=(2,2),
            size = (1000, 500),
            lw=3,
            xlabel="t[s]")
end


# pg = read_powergrid("/Users/jakeroth/git/PowerDynamics.jl/VALIDATION/118bus_fpacopf_SWING.json", Json)
# for k in keys(pg.nodes)
#     node = pg.nodes[k]
#     if typeof(node) == SwingEqLVS
#         pg.nodes[k] = SwingEq(H = node.H, P = node.P, Ω = node.Ω, D = phys[:Dg])
#     end
# end
# write_powergrid(pg, "/Users/jakeroth/git/PowerDynamics.jl/VALIDATION/118bus_fpacopf_SWING.json", Json)
# pg_SWING = read_powergrid("/Users/jakeroth/git/PowerDynamics.jl/VALIDATION/118bus_fpacopf_SWING.json", Json)

# fault96 = LineFailure(line_name=["branch96"],tspan_fault=(1,Inf))
# solution96 = simulate(fault96, pg_SWING, op_pd, (0,1e9));

# solution96 = sol_pd;

# t = solution96.dqsol.t
# u = solution96.dqsol.u
# from_bus = "Bus 38 (load)"
# to_bus = "Bus 65 (generator)"
# gen_bus = "Bus 1 (generator)"
# ft_bus = (38,65)
# g_bus = 1

# cumulative_num_vars = get_cumulative_num_vars(powergrid)
# fbus_idx = cumulative_num_vars[ft_bus[1]]
# tbus_idx = cumulative_num_vars[ft_bus[2]]
# gbus_idx = cumulative_num_vars[g_bus]

# fbus_V = Vector{Float64}(undef, length(t))
# tbus_V = Vector{Float64}(undef, length(t))
# gbus_V = Vector{Float64}(undef, length(t))
# fbus_ang = Vector{Float64}(undef, length(t))
# tbus_ang = Vector{Float64}(undef, length(t))
# gbus_ang = Vector{Float64}(undef, length(t))
# tbus_omega = Vector{Float64}(undef, length(t))
# gbus_omega = Vector{Float64}(undef, length(t))
# for idx in 1:length(t)
#     fbus_V[idx] = abs(u[idx][fbus_idx] + im*u[idx][fbus_idx+1])
#     tbus_V[idx] = abs(u[idx][tbus_idx] + im*u[idx][tbus_idx+1])
#     gbus_V[idx] = abs(u[idx][gbus_idx] + im*u[idx][gbus_idx+1])
#     fbus_ang[idx] = angle(u[idx][fbus_idx] + im*u[idx][fbus_idx+1])
#     tbus_ang[idx] = angle(u[idx][tbus_idx] + im*u[idx][tbus_idx+1])
#     gbus_ang[idx] = angle(u[idx][gbus_idx] + im*u[idx][gbus_idx+1])
#     tbus_omega[idx] = u[idx][tbus_idx+2]
#     gbus_omega[idx] = u[idx][gbus_idx+2]
# end

# v_slice = 1:460
# Plots.plot(t[v_slice], [gbus_V[v_slice], fbus_V[v_slice], tbus_V[v_slice]], 
#     labels = [gen_bus from_bus to_bus],
#     xlabel = "Time (s)",
#     ylabel = "V (p.u.)",
#     seriestype=:line)

# theta_slice = 1:460
# Plots.plot(t[theta_slice], [gbus_ang[theta_slice], fbus_ang[theta_slice], tbus_ang[theta_slice]], 
#     labels = [gen_bus from_bus to_bus],
#     xlabel = "Time (s)",
#     ylabel = "θ (rad)",
#     seriestype=:line)

# omega_slice = 1:460
# Plots.plot(t[omega_slice], [gbus_omega[omega_slice], tbus_omega[omega_slice]], 
#     labels = [gen_bus to_bus],
#     xlabel = "Time (s)",
#     ylabel = "ω (rad/s)",
#     seriestype=:line)