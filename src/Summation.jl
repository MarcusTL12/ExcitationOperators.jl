export summation

function summation(t::CompositeTerm{A}, sum_ind::MOIndex) where {A<:Number}
    deltas = collect(t.deltas)
    is = findall(x -> x.p == sum_ind || x.q == sum_ind, deltas)
    if !isempty(is)
        i = first(is)
        d = deltas[i]
        other_ind = if d.p == sum_ind
            d.q
        else
            d.p
        end

        other_ind_n = other_ind.n

        if other_ind.o == gen
            if sum_ind.o == occ
                other_ind_n = other_ind_n * "ᵒ"
            elseif sum_ind.o == vir
                other_ind_n = other_ind_n * "ᵛ"
            end
        end

        other_ind = ind(sum_ind.o, other_ind_n)

        tensors = copy(t.tensors)
        for i in eachindex(tensors)
            tensors[i] = exchange_index(tensors[i], sum_ind, other_ind)
        end

        operators = copy(t.operators)
        for i in eachindex(operators)
            operators[i] = exchange_index(operators[i], sum_ind, other_ind)
        end

        deltas = SortedSet(
            [d for d in deltas if d.p != sum_ind && d.q != sum_ind]
        )

        CompositeTerm(t.scalar, t.sum_inds, deltas, tensors, operators)
    else
        CompositeTerm(t.scalar, union(t.sum_inds, [sum_ind]),
            t.deltas, t.tensors, t.operators)
    end
end

function summation(s::SumType{A}, sum_ind::MOIndex) where {A<:Number}
    SumType([summation(t, sum_ind) for t in s.terms])
end

function summation(
    a::Union{CompositeTerm{A},SumType{A}},
    sum_inds) where {A<:Number}
    foldl(summation, sum_inds; init=a)
end

