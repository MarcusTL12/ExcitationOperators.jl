export summation

struct SummationType{T<:Number}
    sum_inds::SortedSet{MOIndex}
    term::CompositeTerm{T}
end

function SummationType(sum_inds, s::SummationType{A}) where {A<:Number}
    SummationType(union(sum_inds, s.sum_inds), s.term)
end

function Base.show(io::IO, s::SummationType{T}) where {T<:Number}
    print(io, "Σ_")
    for ind in s.sum_inds
        print(io, ind)
    end

    print(io, '(', s.term, ')')
end

function summation(t::CompositeTerm{A}, sum_ind::MOIndex) where {A<:Number}
    deltas = collect(t.deltas)
    is = findall(x -> x.p == sum_ind || x.q == sum_ind, deltas)
    if !isempty(is)
        i = first(is)
        d = deltas[i]
        other_ind_n = if d.p == sum_ind
            d.q.n
        else
            d.p.n
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

        CompositeTerm(t.scalar, deltas, tensors, operators)
    else
        SummationType(SortedSet(sum_ind), t)
    end
end

function summation(s::SumType{A}, sum_ind::MOIndex) where {A<:Number}
    SumType([summation(t, sum_ind) for t in s.terms])
end

function summation(a::SummationType{A}, sum_ind::MOIndex) where {A<:Number}
    inner_sum = summation(a.term, sum_ind)
    SummationType(a.sum_inds, inner_sum)
end

function summation(
    a::Union{CompositeTerm{A},SumType{A}},
    sum_inds::Vector{MOIndex}) where {A<:Number}
    foldl(summation, sum_inds; init=a)
end

