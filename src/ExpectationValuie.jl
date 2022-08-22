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
                nonop_part = exchange_index(nonop_part, [o.p => np, o.q => nq])
                nonop_part * 2δ(np, nq)
            end
        else
            throw("Have not gotten around to implementing expectation \
            values for strings of $(length(t.operators)) or more excitation \
                operators.")
        end
    end
end
