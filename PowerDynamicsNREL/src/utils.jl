using BlockSystems: namespace_iparams

export BlockPara

"""
    struct BlockPara{P<:Dict}
    BlockPara(block, para; strict=true)

Composite type holdes an `IOBlock` and a parameter `Dict`.
If `strict=true` the construtor asserts that all internal parameters
of the `IOBlock` are given in the parameter dict.
"""
struct BlockPara{P<:Dict}
    block::IOBlock
    para::P
    function BlockPara(block::IOBlock, para::P; strict=true) where {P<:Dict}
        if strict && !(Set(namespace_iparams(block)) ⊆ keys(para))
            throw(ArgumentError("There is a missmatch between given parameter dict and iparams: $(namespace_iparams(block)), $(keys(para))"))
        end
        new{P}(block, para)
    end
end

get_block(bp::BlockPara) = bp.block
get_parameters(bp::BlockPara) = bp.para

"""
    mergep(bps::BlockPara...; strict=true)

Returns tuple of `IOBLock`s ands a merged dict of all parameters.
If `strict=true` check whether the keys in the parameterdicts are unique.
"""
mergep(bps::BlockPara...; kwargs...) = mergep(bps; kwargs...)
function mergep(bps::NTuple{N, BlockPara}; strict=true) where {N}
    isempty(bps) && throw(ArgumentError("Can't merge nothing."))

    paras = get_parameters.(bps)

    if strict && !allunique(vcat(collect.(keys.(paras))...))
        throw(ArgumentError("Keys of parameter dicts not unique!"))
    end

    return get_block.(bps), merge(paras...)
end

subscript(s, i) = Symbol(s, Char(0x02080 + i))
