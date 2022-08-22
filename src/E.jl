struct ExcitationOperator
    p::MOIndex
    q::MOIndex
end

function Base.show(io::IO, e::ExcitationOperator)
    print(io, "E_", e.p.n, e.q.n)
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
