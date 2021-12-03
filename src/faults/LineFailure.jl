
function filter_lines(lines::OrderedDict, line_name)
    OrderedDict(k => v for (k, v) in lines if k != line_name)
end

function filter_lines(lines::Array, line_name)
    copy(lines[1:end.!=line_name])
end

## NEW
function filter_multiple_lines(lines::OrderedDict, line_names)
    OrderedDict(k => v for (k, v) in lines if k ∉ line_names)
end

function filter_multiple_lines(lines::Array, line_names)
    copy(lines[1:end.∉line_names])
end
## NEW

## LineFault -> deprecation

"""
```Julia
LineFault(;from,to)
```
The arguments `from` and `to` specify the line that should be disconnected from the grid.
"""
struct LineFault
    from::Any
    to::Any
    LineFault(; from = from, to = to) = new(from, to)
    @warn "This implementation of a line fault will be deprecated soon. Please use LineFailure instead."
end

function (lf::LineFault)(powergrid)
    @assert lf.from < lf.to "order important to remove the line from powergrid"
    filtered_lines =
        filter(l -> (l.from != lf.from || l.to != lf.to), copy(powergrid.lines))
    PowerGrid(powergrid.nodes, filtered_lines)
end

"""
```Julia
simulate(lf::LineFault, powergrid, x0; timespan)
```
Simulates a [`LineFault`](@ref)
"""
function simulate(lf::LineFault, powergrid, x0; timespan)
    solve(lf(powergrid), x0, timespan)
end

function simulate(lf::LineFault, x0::State; timespan)
    solve(lf(x0.grid), x0.vec, timespan)
end

## NEW
function simulate(lfs::LineFailures, powergrid, x0; timespan)
    solve(lfs(powergrid), x0, timespan)
end

function simulate(lfs::LineFailures, x0::State; timespan)
    solve(lfs(x0.grid), x0.vec, timespan)
end
## NEW

## LineFailure

@doc """
```Julia
LineFailure(;line_name,tspan_fault)
```
The arguments `line_name` and `tspan_fault` specify the line that should be disconnected from the grid and the respective fault duration.
For a list of lines the line_name is the index and for an OrderedDict it is the key of the line.
"""
struct LineFailure <: AbstractPerturbation
    line_name::Any
    tspan_fault::Any
    LineFailure(; line_name = line_name, tspan_fault = tspan_fault) =
        new(line_name, tspan_fault)
end

function (lf::LineFailure)(powergrid)
    filtered_lines = filter_lines(powergrid.lines, lf.line_name)
    PowerGrid(powergrid.nodes, filtered_lines)
end

## NEW
struct LineFailures <: AbstractPerturbation
    line_names::Any
    tspan_fault::Any
    LineFailures(; line_names = line_names, tspan_fault = tspan_fault) =
        new(line_names, tspan_fault)
end

function (lfs::LineFailures)(powergrid)
    filtered_lines = filter_multiple_lines(copy(powergrid.lines), lfs.line_names)
    PowerGrid(powergrid.nodes, filtered_lines)
end
## NEW


export LineFault
export LineFailure
export simulate
