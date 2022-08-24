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

get_indices(t::SymTensor2) = [t.p, t.q]

sym_tensor(symbol, p, q) =
    CompositeTerm(SymTensor2(symbol, p, q))

# Implements the particle-exchange symmetry, g_pqrs = g_rspq
struct SymTensor4 <: Tensor
    symbol::String
    p::MOIndex
    q::MOIndex
    r::MOIndex
    s::MOIndex

    function SymTensor4(symbol, p, q, r, s)
        if (p, q) ≤ (r, s)
            new(symbol, p, q, r, s)
        else
            new(symbol, r, s, p, q)
        end
    end
end

get_indices(t::SymTensor4) = [t.p, t.q, t.r, t.s]

get_symbol(t::T) where {T<:Union{SymTensor2,SymTensor4}} = t.symbol
Base.adjoint(t::T) where {T<:Union{SymTensor2,SymTensor4}} = t

function exchange_index(t::T, i::Int, ind::MOIndex) where
{T<:Union{SymTensor2,SymTensor4}}
    indices = get_indices(t)
    indices[i] = ind
    T(t.symbol, indices...)
end

function exchange_index(t::T, from::MOIndex, to::MOIndex) where
{T<:Union{SymTensor2,SymTensor4}}
    indices = get_indices(t)
    is = findall(x -> x == from, indices)
    for i in is
        indices[i] = to
    end
    T(t.symbol, indices...)
end

sym_tensor(symbol, p, q, r, s) =
    CompositeTerm(SymTensor4(symbol, p, q, r, s))
