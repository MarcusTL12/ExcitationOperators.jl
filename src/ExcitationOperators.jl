module ExcitationOperators

using DataStructures: SortedSet, SortedDict

export occ, vir, gen, ind, E, δ, comm, exval

@enum Occupation gen vir occ

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

function Base.isless(a::ExcitationOperator, b::ExcitationOperator)
    (a.p, a.q) < (b.p, b.q)
end

struct KroeneckerDelta
    p::MOIndex
    q::MOIndex
    function KroeneckerDelta(p::MOIndex, q::MOIndex)
        if p == q
            return 1
        end

        if p > q
            p, q = q, p
        end

        Int(p.o) * Int(q.o) == 2 ? 0 : new(p, q)
    end
end

δ(p::MOIndex, q::MOIndex) = KroeneckerDelta(p, q)

function Base.show(io::IO, δ::KroeneckerDelta)
    print(io, "δ_", δ.p.n, δ.q.n)
end

Base.isless(a::KroeneckerDelta, b::KroeneckerDelta) = (a.p, a.q) < (b.p, b.q)

Base.adjoint(d::KroeneckerDelta) = d

struct OperatorProduct
    deltas::SortedSet{KroeneckerDelta}
    operators::Vector{ExcitationOperator}
end

function Base.:(==)(a::OperatorProduct, b::OperatorProduct)
    (collect(a.deltas), a.operators) == (collect(b.deltas), b.operators)
end

function Base.:*(a::KroeneckerDelta, b::KroeneckerDelta)
    OperatorProduct(SortedSet([a, b]), ExcitationOperator[])
end

function Base.:*(a::OperatorProduct, b::KroeneckerDelta)
    OperatorProduct(union(a.deltas, [b]), a.operators)
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

function Base.isless(a::OperatorProduct, b::OperatorProduct)
    (a.operators, a.deltas) < (b.operators, b.deltas)
end

function Base.isless(a::ExcitationOperator, b::OperatorProduct)
    ([a], KroeneckerDelta[]) < (b.operators, b.deltas)
end
function Base.isless(a::OperatorProduct, b::ExcitationOperator)
    (a.operators, a.deltas) < ([b], KroeneckerDelta[])
end

function Base.isless(a::KroeneckerDelta, b::OperatorProduct)
    (ExcitationOperator[], [a]) < (b.operators, b.deltas)
end
function Base.isless(a::OperatorProduct, b::KroeneckerDelta)
    (a.operators, a.deltas) < (ExcitationOperator[], [b])
end

Base.isless(::KroeneckerDelta, ::ExcitationOperator) = true
Base.isless(::ExcitationOperator, ::KroeneckerDelta) = false

Base.isless(::Nothing, ::KroeneckerDelta) = true
Base.isless(::KroeneckerDelta, ::Nothing) = false

Base.isless(::Nothing, ::OperatorProduct) = true
Base.isless(::OperatorProduct, ::Nothing) = false

Base.isless(::Nothing, ::ExcitationOperator) = true
Base.isless(::ExcitationOperator, ::Nothing) = false

Base.isless(::Nothing, ::Nothing) = false

const OpUnion =
    Union{Nothing,KroeneckerDelta,ExcitationOperator,OperatorProduct}

struct OperatorSum{T<:Number}
    s::Vector{Pair{OpUnion,T}}

    function OperatorSum(s::Vector{Pair{OpUnion,T}}) where {T<:Number}
        lookup = SortedDict{OpUnion,T}()
        for (o, n) in s
            lookup[o] = get(lookup, o, zero(T)) + n
        end
        s2 = Pair{OpUnion,T}[o => n for (o, n) in lookup if !iszero(n)]
        isempty(s2) ? zero(T) : new{T}(s2)
    end
end

function Base.show(io::IO, s::OperatorSum)
    if !isone(abs(s.s[1].second))
        print(io,
            s.s[1].second, ' ', isnothing(s.s[1].first) ? '\b' : s.s[1].first)
    else
        if isone(-s.s[1].second)
            print(io, '-')
        end
        print(io, isnothing(s.s[1].first) ? 1 : s.s[1].first)
    end
    for o in @view s.s[2:end]
        print(io, ' ', o.second < 0 ? '-' : '+', ' ')
        if !isone(abs(o.second))
            print(io, abs(o.second), ' ')
        end
        print(io, isnothing(o.first) ? 1 : o.first)
    end
end

function Base.:-(a::A) where {A<:OpUnion}
    OperatorSum(Pair{OpUnion,Int}[a=>-1])
end

function Base.:*(n::T, p::OP) where {T<:Number,OP<:OpUnion}
    OperatorSum(Pair{OpUnion,T}[p=>n])
end
Base.:*(p::OP, n::T) where {T<:Number,OP<:OpUnion} = n * p

function Base.:+(a::A, b::B) where {A<:OpUnion,B<:OpUnion}
    OperatorSum(Pair{OpUnion,Int}[a=>1, b=>1])
end
function Base.:-(a::A, b::B) where {A<:OpUnion,B<:OpUnion}
    OperatorSum(Pair{OpUnion,Int}[a=>1, b=>-1])
end

function Base.:+(a::A, b::B) where {A<:Number,B<:OpUnion}
    OperatorSum(Pair{OpUnion,A}[nothing=>a, b=>one(A)])
end
function Base.:+(a::A, b::B) where {A<:OpUnion,B<:Number}
    b + a
end

Base.:*(::Nothing, b::B) where {B<:OpUnion} = b
Base.:*(a::A, ::Nothing) where {A<:OpUnion} = a

function Base.:-(a::A, b::B) where {A<:Number,B<:OpUnion}
    OperatorSum(Pair{OpUnion,A}[nothing=>a, b=>-one(A)])
end
function Base.:-(a::A, b::B) where {A<:OpUnion,B<:Number}
    b + (-a)
end

function Base.:+(a::OperatorSum{T}, b::OP) where {T<:Number,OP<:OpUnion}
    OperatorSum(Pair{OpUnion,Int}[a.s; b => 1])
end
function Base.:+(a::OP, b::OperatorSum{T}) where {T<:Number,OP<:OpUnion}
    OperatorSum(Pair{OpUnion,promote_type(Int, T)}[a => 1; b.s])
end

function Base.:+(
    a::OperatorSum{A}, b::OperatorSum{B}
) where {A<:Number,B<:Number}
    OperatorSum(Pair{OpUnion,promote_type(A, B)}[a.s; b.s])
end

function Base.:+(a::A, b::OperatorSum{B}) where {A<:Number,B<:Number}
    OperatorSum(Pair{OpUnion,promote_type(A, B)}[nothing => a; b.s])
end

function Base.:+(a::OperatorSum{A}, b::B) where {A<:Number,B<:Number}
    OperatorSum(Pair{OpUnion,promote_type(A, B)}[a.s; nothing => b])
end

function Base.:-(a::OperatorSum{T}) where {T<:Number}
    OperatorSum(Pair{OpUnion,T}[o => -n for (o, n) in a.s])
end

function Base.:-(a::OperatorSum{T}, b::OP) where {T<:Number,OP<:OpUnion}
    OperatorSum(Pair{OpUnion,promote_type(Int, T)}[a.s; b => -1])
end

function Base.:-(a::OP, b::OperatorSum{T}) where {T<:Number,OP<:OpUnion}
    a + (-b)
end

function Base.:-(
    a::OperatorSum{A}, b::OperatorSum{B}
) where {A<:Number,B<:Number}
    a + (-b)
end

function Base.:-(a::A, b::OperatorSum{B}) where {A<:Number,B<:Number}
    a + (-b)
end
function Base.:-(a::OperatorSum{A}, b::B) where {A<:Number,B<:Number}
    a + (-b)
end

function Base.:*(a::A, b::OperatorSum{B}) where {A<:Number,B<:Number}
    OperatorSum(Pair{OpUnion,promote_type(A, B)}[o => n * a for (o, n) in b.s])
end
function Base.:*(a::OperatorSum{A}, b::B) where {A<:Number,B<:Number}
    b * a
end

Base.adjoint(::Nothing) = nothing

function Base.adjoint(a::OperatorSum{T}) where {T<:Number}
    OperatorSum(Pair{OpUnion,T}[o' => n' for (o, n) in a.s])
end

function Base.:*(
    a::OperatorSum{A}, b::OperatorSum{B}
) where {A<:Number,B<:Number}
    OperatorSum(
        Pair{OpUnion,promote_type(A, B)}[
            o1 * o2 => n1 * n2
            for (o1, n1) in a.s, (o2, n2) in b.s
        ][:]
    )
end

function Base.:*(a::OperatorSum{A}, b::B) where {A<:Number,B<:OpUnion}
    OperatorSum(Pair{OpUnion,A}[o * b => n for (o, n) in a.s])
end
function Base.:*(a::A, b::OperatorSum{B}) where {A<:OpUnion,B<:Number}
    OperatorSum(Pair{OpUnion,B}[a * o => n for (o, n) in b.s])
end

comm(a, b::T) where {T<:Number} = 0
comm(a::T, b) where {T<:Number} = 0

function comm(a::ExcitationOperator, b::ExcitationOperator)
    δ(a.q, b.p) * E(a.p, b.q) - δ(a.p, b.q) * E(b.p, a.q)
end

Base.one(::Type{ExcitationOperator}) = 1
Base.one(::Type{KroeneckerDelta}) = 1

function comm(a::ExcitationOperator, b::OperatorProduct)
    prod(b.deltas) * sum(
        prod(b.operators[1:i-1]) *
        comm(a, b.operators[i]) *
        prod(b.operators[i+1:end])
        for i in 1:length(b.operators)
    )
end

function comm(
    a::OperatorProduct, b::B
) where {B<:Union{ExcitationOperator,OperatorProduct}}
    prod(a.deltas) * sum(
        prod(a.operators[1:i-1]) *
        comm(a.operators[i], b) *
        prod(a.operators[i+1:end])
        for i in 1:length(a.operators)
    )
end

function comm(
    a::A, b::OperatorSum{B}
) where {A<:Union{ExcitationOperator,OperatorProduct},B<:Number}
    sum(n * comm(a, o) for (o, n) in b.s)
end

function comm(a::OperatorSum{A}, b) where {A<:Number}
    sum(n * comm(o, b) for (o, n) in a.s)
end

# Basic commutator arithmetic works. Here we start expectation values.

struct ExpectationValue
    p::OperatorProduct
end

function Base.:(==)(a::ExpectationValue, b::ExpectationValue)
    a.p == b.p
end

function Base.show(io::IO, ev::ExpectationValue)
    printed = false

    for d in ev.p.deltas
        printed = true
        print(io, d, ' ')
    end

    print(io, '⟨')

    for e in ev.p.operators
        print(io, e, ' ')
    end

    print(io, '\b')

    print(io, '⟩')
end

exval(n::T) where {T<:Number} = n
exval(::Nothing) = nothing
exval(d::KroeneckerDelta) = d

function exval(e::ExcitationOperator)
    if e.p.o == vir || e.q.o == vir
        0
    elseif e.p.o == occ || e.q.o == occ
        2δ(e.p, e.q)
    else
        ExpectationValue(OperatorProduct(SortedSet(), [e]))
    end
end

function exval(p::OperatorProduct)
    if isempty(p.operators)
        p
    elseif p.operators[1].p.o == vir || p.operators[end].q.o == vir
        0
    elseif p.operators[1].q.o == vir
        exval(comm(
            p.operators[1], OperatorProduct(p.deltas, p.operators[2:end])
        ))
    elseif p.operators[end].p.o == vir
        exval(comm(
            OperatorProduct(p.deltas, p.operators[1:end-1]), p.operators[end]
        ))
    else
        ExpectationValue(p)
    end
end

Base.isless(a::ExpectationValue, b::ExpectationValue) = a.p < b.p

struct ExpectationProduct
    vals::Vector{ExpectationValue}
    function ExpectationProduct(v::Vector{ExpectationValue})
        v = sort(v)
        for i = 2:length(v)
            union!(v[1].p.deltas, v[i].p.deltas)
            empty!(v[i].p.deltas)
        end
        new(v)
    end
end

function Base.show(io::IO, p::ExpectationProduct)
    for v in p.vals
        print(io, v)
    end
end

function Base.:*(a::ExpectationValue, b::ExpectationValue)
    ExpectationProduct([a, b])
end

const ExUnion = 0

struct ExpectationSum{T<:Number}

end

end # module
