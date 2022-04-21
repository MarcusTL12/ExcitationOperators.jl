module ExcitationOperators

using DataStructures: SortedSet

export Occupation, MOIndex, ExcitationOperator, KroeneckerDelta, occ, vir, gen

export ind, E, δ

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

end # module
