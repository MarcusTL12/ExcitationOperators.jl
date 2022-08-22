export exval

function exval(t::CompositeTerm{T}) where {T<:Number}
    if isempty(t.operators)
        t
    else
        if first(t.operators).p.o == vir || last(t.operators).q.o == vir
            return CompositeTerm(zero(T))
        end

        nonop_part = get_nonop(t)

        if length(t.operators) == 1
            o = first(t.operators)
            if o.p.o == occ || o.q.o == occ
                nonop_part * 2δ(o.p, o.q)
            else
                np = make_occ(o.p)
                nq = make_occ(o.q)
                exchange_index(nonop_part, [o.p => np, o.q => nq]) * 2δ(np, nq)
            end
        else
            i = findfirst(x -> x.p.o == vir || x.q.o == vir, t.operators)
            if i isa Int
                if t.operators[i].p.o == vir
                    i -= 1
                end

                 exval(nonop_part * comm(
                    CompositeTerm(t.operators[1:i]),
                    CompositeTerm(t.operators[i+1:end])
                ))
            elseif first(t.operators).q.o == occ
                nonop_part * exval(CompositeTerm(first(t.operators))) *
                exval(CompositeTerm(t.operators[2:end]))
            elseif last(t.operators).p.o == occ
                nonop_part * exval(CompositeTerm(t.operators[1:end-1])) *
                exval(CompositeTerm(last(t.operators)))
            else
                throw("Not implemented: $t")
            end
        end
    end
end

function exval(s::SumType{T}) where {T<:Number}
    SumType(exval.(s.terms))
end
