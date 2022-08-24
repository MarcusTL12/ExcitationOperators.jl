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

        p.o == occ && q.o == vir || p.o == vir && q.o == occ ? 0 : new(p, q)
    end
end

function Base.show(io::IO, δ::KroeneckerDelta)
    print(io, "δ_", δ.p.n, δ.q.n)
end

Base.isless(a::KroeneckerDelta, b::KroeneckerDelta) = (a.p, a.q) < (b.p, b.q)

Base.adjoint(d::KroeneckerDelta) = d

function exchange_index(d::KroeneckerDelta, i::Int, ind::MOIndex)
    if i == 1
        KroeneckerDelta(ind, d.q)
    elseif i == 2
        KroeneckerDelta(d.p, ind)
    end
end

function exchange_index(d::KroeneckerDelta, from::MOIndex, to::MOIndex)
    KroeneckerDelta(d.p == from ? to : d.p, d.q == from ? to : d.q)
end

export δ
δ(p::MOIndex, q::MOIndex) = CompositeTerm(KroeneckerDelta(p, q))
