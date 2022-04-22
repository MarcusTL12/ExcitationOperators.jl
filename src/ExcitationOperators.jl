module ExcitationOperators

using DataStructures: SortedSet

export Occupation, MOIndex, ExcitationOperator, KroeneckerDelta, occ, vir, gen

export ind, E, δ, comm

@enum Occupation gen occ vir

struct MOIndex
    o::Occupation
    n::String
end

ind(o, name) = MOIndex(o, name)

Base.isless(p::MOIndex, q::MOIndex) = (p.o, p.n) < (q.o, q.n)

struct ExcitationOperator
    p::MOIndex
    q::MOIndex
end

E(p::MOIndex, q::MOIndex) = ExcitationOperator(p, q)

function Base.show(io::IO, e::ExcitationOperator)
    print(io, "E_", e.p.n, e.q.n)
end

Base.adjoint(e::ExcitationOperator) = E(e.q, e.p)

struct KroeneckerDelta
    p::MOIndex
    q::MOIndex
end

function δ(p::MOIndex, q::MOIndex)
    if p > q
        p, q = q, p
    end
    KroeneckerDelta(p, q)
end

function Base.show(io::IO, δ::KroeneckerDelta)
    print(io, "δ_", δ.p.n, δ.q.n)
end

Base.isless(a::KroeneckerDelta, b::KroeneckerDelta) = (a.p, a.q) < (b.p, b.q)

struct OperatorProduct
    deltas::SortedSet{KroeneckerDelta}
    operators::Vector{ExcitationOperator}
end

function Base.:*(a::KroeneckerDelta, b::KroeneckerDelta)
    OperatorProduct(SortedSet([a, b]), ExcitationOperator[])
end

function Base.:*(a::OperatorProduct, b::KroeneckerDelta)
    OperatorProduct(union(a.deltas, b), a.operators)
end

Base.:*(a::KroeneckerDelta, b::OperatorProduct) = b * a

function Base.:*(a::ExcitationOperator, b::ExcitationOperator)
    OperatorProduct(SortedSet(), ExcitationOperator[a, b])
end

function Base.:*(a::KroeneckerDelta, b::ExcitationOperator)
    OperatorProduct(SortedSet([a]), ExcitationOperator[b])
end

Base.:*(a::ExcitationOperator, b::KroeneckerDelta) = b * a

function Base.:*(a::OperatorProduct, b::ExcitationOperator)
    OperatorProduct(a.deltas, vcat(a.operators, b))
end

function Base.:*(a::ExcitationOperator, b::OperatorProduct)
    OperatorProduct(b.deltas, vcat(a, b.operators))
end

function Base.:*(a::OperatorProduct, b::OperatorProduct)
    OperatorProduct(union(a.deltas, b.deltas), vcat(a.operators, b.operators))
end

function Base.adjoint(p::OperatorProduct)
    OperatorProduct(p.deltas, reverse!([x' for x in p.operators]))
end

function Base.show(io::IO, p::OperatorProduct)
    printed = false

    for d in p.deltas
        printed = true
        print(io, d, ' ')
    end

    for e in p.operators
        printed = true
        print(io, e, ' ')
    end

    if printed
        print(io, '\b')
    end
end

const OpUnion =
    Union{Nothing,KroeneckerDelta,ExcitationOperator,OperatorProduct}

struct OperatorSum{T<:Number}
    s::Vector{Pair{T,OpUnion}}
end

function Base.show(io::IO, s::OperatorSum)
    if !isone(abs(s.s[1].first))
        print(io, s.s[1].first, ' ')
    elseif isone(-s.s[1].first)
        print(io, '-')
    end
    print(io, isnothing(s.s[1].second) ? '\b' : s.s[1].second)
    for o in @view s.s[2:end]
        print(io, ' ', o.first < 0 ? '-' : '+', ' ')
        if !isone(abs(o.first))
            print(io, abs(o.first), ' ')
        end
        print(io, isnothing(o.second) ? '\b' : o.second)
    end
end

function Base.:*(n::T, p::OP) where {T<:Number,OP<:OpUnion}
    OperatorSum(Pair{T,OpUnion}[n=>p])
end

Base.:*(p::OP, n::T) where {T<:Number,OP<:OpUnion} = n * p

function Base.:+(a::A, b::B) where {A<:OpUnion,B<:OpUnion}
    OperatorSum(Pair{Int,OpUnion}[1=>a, 1=>b])
end

function Base.:-(a::A, b::B) where {A<:OpUnion,B<:OpUnion}
    OperatorSum(Pair{Int,OpUnion}[1=>a, -1=>b])
end

function Base.:+(a::OperatorSum{T}, b::OP) where {T<:Number,OP<:OpUnion}
    OperatorSum(Pair{Int,OpUnion}[a.s; 1 => b])
end

function Base.:+(a::OP, b::OperatorSum{T}) where {T<:Number,OP<:OpUnion}
    OperatorSum(Pair{promote_type(Int, T),OpUnion}[1 => a; b.s])
end

function Base.:+(
    a::OperatorSum{A}, b::OperatorSum{B}
) where {A<:Number,B<:Number}
    OperatorSum(Pair{promote_type(A, B),OpUnion}[a.s; b.s])
end

function Base.:-(a::OperatorSum{T}) where {T<:Number}
    OperatorSum(Pair{T,OpUnion}[-n => o for (n, o) in a.s])
end

function Base.:-(a::OperatorSum{T}, b::OP) where {T<:Number,OP<:OpUnion}
    OperatorSum(Pair{promote_type(Int, T),OpUnion}[a.s; -1 => b])
end

function Base.:-(a::OP, b::OperatorSum{T}) where {T<:Number,OP<:OpUnion}
    a + (-b)
end

function Base.:-(
    a::OperatorSum{A}, b::OperatorSum{B}
) where {A<:Number,B<:Number}
    a + (-b)
end

function Base.:*(a::A, b::OperatorSum{B}) where {A<:Number,B<:Number}
    OperatorSum(Pair{promote_type(A, B),OpUnion}[n * a => o for (n, o) in b.s])
end

function Base.:*(a::OperatorSum{A}, b::B) where {A<:Number,B<:Number}
    b * a
end

# function comm(a::ExcitationOperator, b::ExcitationOperator)
#     d(a.q, b.q) * E(a.p, b.q) - d(a.p, b.q) * E(b.p, a.q)
# end

end # module
