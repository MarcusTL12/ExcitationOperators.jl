export summation

function summation(t::CompositeTerm{A}, sum_ind::MOIndex) where {A<:Number}
    deltas = collect(t.deltas)
    is = findall(x -> x.p == sum_ind || x.q == sum_ind, deltas)
    if !isempty(is)
        i = first(is)
        d = deltas[i]
        old_ind = if d.p == sum_ind
            d.q
        else
            d.p
        end

        new_ind_n = old_ind.n

        # if old_ind.o == gen
        #     if sum_ind.o == occ
        #         new_ind_n = new_ind_n
        #     elseif sum_ind.o == vir
        #         new_ind_n = new_ind_n
        #     end
        # end

        new_ind_o = old_ind.o

        if sum_ind.o != gen
            new_ind_o = sum_ind.o
        end

        new_ind = ind(new_ind_o, new_ind_n)

        tensors = copy(t.tensors)
        for i in eachindex(tensors)
            tensors[i] = exchange_index(tensors[i], sum_ind, new_ind)
            tensors[i] = exchange_index(tensors[i], old_ind, new_ind)
        end

        operators = copy(t.operators)
        for i in eachindex(operators)
            operators[i] = exchange_index(operators[i], sum_ind, new_ind)
            operators[i] = exchange_index(operators[i], old_ind, new_ind)
        end

        new_deltas = SortedSet{KroeneckerDelta}()
        for d in deltas
            nd = exchange_index(d, sum_ind, new_ind)
            if nd isa KroeneckerDelta
                nd = exchange_index(nd, old_ind, new_ind)
            end
            if nd isa KroeneckerDelta
                push!(new_deltas, nd)
            elseif iszero(nd)
                return CompositeTerm(zero(A))
            end
        end

        sum_inds = collect(t.sum_inds)
        for i in eachindex(sum_inds)
            if sum_inds[i] == sum_ind || sum_inds[i] == old_ind
                sum_inds[i] == new_ind
            end
        end
        sum_inds = SortedSet(sum_inds)

        CompositeTerm(t.scalar, sum_inds, new_deltas, tensors, operators)
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

