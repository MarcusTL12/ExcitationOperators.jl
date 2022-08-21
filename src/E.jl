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
