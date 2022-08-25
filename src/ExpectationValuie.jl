export exval

function exval(t::CompositeTerm{T}) where {T<:Number}
    if isempty(t.operators)
        t
    else
        if first(t.operators).p.o == vir || last(t.operators).q.o == vir
            return CompositeTerm(zero(T))
        end

        first_ind = first(t.operators).p
        last_ind = last(t.operators).q

        t = exchange_index(
            t,
            [
                first_ind => make_occ(first_ind),
                last_ind => make_occ(last_ind)
            ]
        )

        nonop_part = get_nonop(t)

        if length(t.operators) == 1
            o = first(t.operators)
            if o.p.o == occ || o.q.o == occ
                summation(get_nosum(nonop_part) * 2δ(o.p, o.q),
                    nonop_part.sum_inds)
            else
                np = make_occ(o.p)
                nq = make_occ(o.q)
                nonop_part = exchange_index(nonop_part, [o.p => np, o.q => nq])
                summation(get_nosum(nonop_part) * 2δ(np, nq),
                    nonop_part.sum_inds)
            end
        else
            i = findfirst(x -> x.p.o == vir || x.q.o == vir, t.operators)
            if i isa Int
                if t.operators[i].p.o == vir
                    i -= 1
                end

                exval(summation(get_nosum(nonop_part) * comm(
                        CompositeTerm(t.operators[1:i]),
                        CompositeTerm(t.operators[i+1:end])
                    ), nonop_part.sum_inds))
            elseif first(t.operators).q.o == occ
                exval(CompositeTerm(first(t.operators))) *
                exval(summation(
                    get_nosum(nonop_part) * CompositeTerm(t.operators[2:end]),
                    nonop_part.sum_inds
                ))
            elseif last(t.operators).p.o == occ
                exval(summation(
                    get_nosum(nonop_part) * CompositeTerm(t.operators[1:end-1]),
                    nonop_part.sum_inds
                )) *
                exval(CompositeTerm(last(t.operators)))
            else
                o = first(t.operators)
                q_occ = make_occ(o.q)
                q_vir = make_vir(o.q)

                t_occ = exchange_index(t, o.q, q_occ)
                t_vir = exchange_index(t, o.q, q_vir)

                # @show t_occ t_vir

                exval(t_occ) + exval(t_vir)
            end
        end
    end
end

function exval(s::SumType{T}) where {T<:Number}
    sum(exval(t) for t in s.terms)
end
