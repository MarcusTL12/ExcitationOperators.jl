# File for hermitian style symmetry: hpq <-> hqp, gpqrs <-> gqpsr
export sym_tensor

# Implements symmetric 2-tensor h_pq = h_qp
struct SymTensor2 <: Tensor
    symbol::String
    p::MOIndex
    q::MOIndex

    function SymTensor2(symbol, p, q)
        if p <= q
            new(symbol, p, q)
        else
            new(symbol, q, p)
        end
    end
end

get_symbol(t::T) where {T<:SymTensor2} = t.symbol
get_indices(t::SymTensor2) = [t.p, t.q]
Base.adjoint(t::T) where {T<:SymTensor2} = t

function exchange_index(t::SymTensor2, i::Int, ind::MOIndex)
    indices = get_indices(t)
    indices[i] = ind
    SymTensor2(t.symbol, indices...)
end

function exchange_index(t::SymTensor2, from::MOIndex, to::MOIndex)
    indices = get_indices(t)
    is = findall(x -> x == from, indices)
    for i in is
        indices[i] = to
    end
    SymTensor2(t.symbol, indices...)
end

sym_tensor(symbol, p, q) =
    CompositeTerm(SymTensor2(symbol, p, q))
