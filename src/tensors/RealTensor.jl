export real_tensor

struct RealTensor <: Tensor
    symbol::String
    indices::Vector{MOIndex}
end

get_symbol(t::RealTensor) = t.symbol
get_indices(t::RealTensor) = t.indices
Base.adjoint(t::RealTensor) = t

function exchange_index(t::RealTensor, i::Int, ind::MOIndex)
    indices = copy(t.indices)
    indices[i] = ind
    RealTensor(t.symbol, indices)
end

function exchange_index(t::RealTensor, from::MOIndex, to::MOIndex)
    indices = copy(t.indices)
    is = findall(x -> x == from, indices)
    for i in is
        indices[i] = to
    end
    RealTensor(t.symbol, indices)
end

real_tensor(symbol, indices...) =
    CompositeTerm(RealTensor(symbol, collect(indices)))
