# CC amplitude tensor printed as t_ijkl^abcd and appropriate latex code
# wrapper around PSymTensor, so implements particle exchange symmetry by default
# cc_amp_tensor("t", a, i, b, j) -> t_ij^ab
export cc_amp_tensor

struct CCAmpTensor <: Tensor
    inner::PSymTensor

    function CCAmpTensor(symbol::String, indices::Vector{MOIndex})
        new(PSymTensor(symbol, indices))
    end
end

get_symbol(t::CCAmpTensor) = t.inner.symbol
get_indices(t::CCAmpTensor) = t.inner.indices
Base.adjoint(t::CCAmpTensor) = t

function exchange_index(t::CCAmpTensor, i::Int, ind::MOIndex)
    indices = copy(t.inner.indices)
    indices[i] = ind
    CCAmpTensor(t.inner.symbol, indices)
end

function exchange_index(t::CCAmpTensor, from::MOIndex, to::MOIndex)
    indices = copy(t.inner.indices)
    is = findall(x -> x == from, indices)
    for i in is
        indices[i] = to
    end
    CCAmpTensor(t.inner.symbol, indices)
end

cc_amp_tensor(symbol, indices...) =
    CompositeTerm(CCAmpTensor(symbol, collect(indices)))

function Base.show(io::IO, t::CCAmpTensor)
    print(io, get_symbol(t), '_')
    inds = get_indices(t)
    for ind in inds[2:2:end]
        print(io, ind)
    end
    print(io, '^')
    for ind in inds[1:2:end]
        print(io, ind)
    end
end

function printlatex(io::IO, t::CCAmpTensor, color)
    print(io, get_symbol(t), "_{")
    inds = get_indices(t)
    printlatex(io, inds[2], color)
    for i in inds[4:2:end]
        print(io, ' ')
        printlatex(io, i, color)
    end
    print(io, "}^{")
    printlatex(io, inds[1], color)
    for i in inds[3:2:end]
        print(io, ' ')
        printlatex(io, i, color)
    end
    print(io, '}')
end
