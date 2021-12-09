function simple_energy_model_failure(opfdata, case_options, optimal_values, l)
    opfd = deepcopy(opfdata)
    nline = length(opfd.lines)
    opfd.lines = opfd.lines[filter(x->x!=l,1:nline)]
    opfmodeldata = get_opfmodeldata(opfd, case_options);
    opfmodeldata[:Y] = imag.(opfmodeldata[:Y])
    op, status = simple_energy_model(opfmodeldata, optimal_values)
    return op, status
end

function simple_energy_model(opfmodeldata, optimal_values)
    buses        = opfmodeldata[:buses]
    nonLoadBuses = opfmodeldata[:nonLoadBuses]
    bus_ref      = opfmodeldata[:bus_ref]
    Y            = opfmodeldata[:Y]
    nbus         = length(buses)
    nrow         = 2nbus - length(nonLoadBuses) - 1
    VMbar        = optimal_values[:Vm]
    VAbar        = optimal_values[:Va]
    baseMVA      = opfmodeldata[:baseMVA]

    em = Model(Ipopt.Optimizer)

    #
    # Variables
    #
    @variable(em, Vm[1:nbus])
    @variable(em, -pi <= Va[1:nbus] <= pi)
    setlowerbound(Va[bus_ref], VAbar[bus_ref])
    setupperbound(Va[bus_ref], VAbar[bus_ref])
    setlowerbound(Vm[bus_ref], VMbar[bus_ref])
    setupperbound(Vm[bus_ref], VMbar[bus_ref])
    for b in 1:nbus
        if b in nonLoadBuses
            setlowerbound(Vm[b], VMbar[b])
            setupperbound(Vm[b], VMbar[b])
        end
    end

    #
    # Set objective to minimize energy function
    #
    @NLobjective(em, Min,
        -0.5*sum(Y[m,n]*Vm[m]*Vm[n]*cos(Va[m] - Va[n]) for m in 1:nbus for n in 1:nbus if Y[m,n] != 0)
        - sum(optimal_values[:Pnet][m]*Va[m] + optimal_values[:Qnet][m]*log(Vm[m]) for m in 1:nbus))

    #
    # initial value
    #
    setvalue.(getindex(em, :Vm), VMbar)
    setvalue.(getindex(em, :Va), VAbar)


    #
    # solve model
    #
    optimize!(em)
    out = Dict()
    out[:Vm] = value.(Vm)
    out[:Va] = value.(Va)
    status = termination_status(em)
    return (out, status)
end