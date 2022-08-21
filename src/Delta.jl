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

function Base.show(io::IO, δ::KroeneckerDelta)
    print(io, "δ_", δ.p.n, δ.q.n)
end

Base.isless(a::KroeneckerDelta, b::KroeneckerDelta) = (a.p, a.q) < (b.p, b.q)

Base.adjoint(d::KroeneckerDelta) = d
