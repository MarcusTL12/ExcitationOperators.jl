abstract type Tensor end

# These methods should be implemented for all subtypes of Tensor

get_symbol(::T) where {T<:Tensor} =
    throw("get_symbol not implemented for Tensor type $(T)!")
get_indices(::T) where {T<:Tensor} =
    throw("get_indices not implemented for Tensor type $(T)!")
Base.adjoint(::T) where {T<:Tensor} =
    throw("adjoint not implemented for Tensor type $(T)!")
exchange_index(::T, i::Int, ind::MOIndex) where {T<:Tensor} =
    throw("exchange_index not implemented for Tensor type $(T)")
exchange_index(::T, from::MOIndex, to::MOIndex) where {T<:Tensor} =
    throw("exchange_index not implemented for Tensor type $(T)")

# Base.show is overridable if wanted (typically for expectation values)
function Base.show(io::IO, t::T) where {T<:Tensor}
    print(io, get_symbol(t), '_')
    for ind in get_indices(t)
        print(io, ind.n)
    end
end

function Base.:(==)(a::A, b::B) where {A<:Tensor,B<:Tensor}
    (get_symbol(a), get_indices(a)) == (get_symbol(b), get_indices(b))
end

function Base.isless(a::A, b::B) where {A<:Tensor,B<:Tensor}
    (get_symbol(a), get_indices(a)) < (get_symbol(b), get_indices(b))
end

function Base.hash(t::T) where {T<:Tensor}
    hash((get_symbol(t), get_indices(t)))
end

function Base.hash(v::Vector{Tensor})
    if isempty(v)
        hash(314)
    else
        h = hash(first(v))
        for t in v[2:end]
            h = hash(h * hash(t))
        end
        h
    end
end

# Basic tensor implementation:

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


struct SymTensor4 <: Tensor
    # Implements the particle-exchange symmetry, g_pqrs = g_rspq
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

get_symbol(t::SymTensor4) = t.symbol
get_indices(t::SymTensor4) = [t.p, t.q, t.r, t.s]
Base.adjoint(t::SymTensor4) = t

function exchange_index(t::SymTensor4, i::Int, ind::MOIndex)
    indices = get_indices(t)
    indices[i] = ind
    SymTensor4(t.symbol, indices...)
end

function exchange_index(t::SymTensor4, from::MOIndex, to::MOIndex)
    indices = get_indices(t)
    is = findall(x -> x == from, indices)
    for i in is
        indices[i] = to
    end
    SymTensor4(t.symbol, indices...)
end
