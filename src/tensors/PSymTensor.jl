# File for particle exchange symmetry: g(pq)(rs) <-> g(rs)(pq)
export psym_tensor

struct PSymTensor <: Tensor
    symbol::String
    indices::Vector{MOIndex}

    function PSymTensor(symbol::String, indices::Vector{MOIndex})
        if !iseven(length(indices))
            throw(
                "Particle exchange symmetric tensors must have even dimension"
            )
        end

        new(
            symbol,
            reshape(
                sortslices(reshape(indices, 2, length(indices) ÷ 2), dims=2),
                length(indices)
            )
        )
    end
end

get_symbol(t::PSymTensor) = t.symbol
get_indices(t::PSymTensor) = t.indices
Base.adjoint(t::PSymTensor) = t

function exchange_index(t::PSymTensor, i::Int, ind::MOIndex)
    indices = copy(t.indices)
    indices[i] = ind
    PSymTensor(t.symbol, indices)
end

function exchange_index(t::PSymTensor, from::MOIndex, to::MOIndex)
    indices = copy(t.indices)
    is = findall(x -> x == from, indices)
    for i in is
        indices[i] = to
    end
    PSymTensor(t.symbol, indices)
end

psym_tensor(symbol, indices...) =
    CompositeTerm(PSymTensor(symbol, collect(indices)))

