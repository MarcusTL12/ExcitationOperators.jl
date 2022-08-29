export E, e

struct ExcitationOperator
    p::MOIndex
    q::MOIndex
end

function Base.show(io::IO, e::ExcitationOperator)
    print(io, "E_", e.p, e.q)
end

Base.adjoint(e::ExcitationOperator) = ExcitationOperator(e.q, e.p)

function Base.isless(a::ExcitationOperator, b::ExcitationOperator)
    (a.p, a.q) < (b.p, b.q)
end

function exchange_index(o::ExcitationOperator, i::Int, ind::MOIndex)
    if i == 1
        ExcitationOperator(ind, o.q)
    elseif i == 2
        ExcitationOperator(o.p, ind)
    end
end

function exchange_index(o::ExcitationOperator, from::MOIndex, to::MOIndex)
    ExcitationOperator(o.p == from ? to : o.p, o.q == from ? to : o.q)
end

E(p::MOIndex, q::MOIndex) = CompositeTerm(ExcitationOperator(p, q))

e(p::MOIndex, q::MOIndex, r::MOIndex, s::MOIndex) =
    E(p, q) * E(r, s) - δ(q, r) * E(p, s)
