export comm

function comm(::A, ::B) where {A<:Number,B}
    zero(promote_type(A, B))
end

function comm(::A, ::B) where {A,B<:Number}
    zero(promote_type(A, B))
end

function comm(a::ExcitationOperator, b::ExcitationOperator)
    δ(a.q, b.p) * E(a.p, b.q) - δ(a.p, b.q) * E(b.p, a.q)
end

function comm(a::ExcitationOperator, b::Vector{ExcitationOperator})
    if length(b) == 1
        comm(a, b[1])
    elseif length(b) == 2
        CompositeTerm(b[1]) * comm(a, b[2]) +
        comm(a, b[1]) * CompositeTerm(b[2])
    else
        sum(
            prod(CompositeTerm, b[1:i-1]; init=CompositeTerm(1)) *
            comm(a, b[i]) *
            prod(CompositeTerm, b[i+1:end]; init=CompositeTerm(1))
            for i in eachindex(b)
        )
    end
end

function comm(a::Vector{ExcitationOperator}, b::Vector{ExcitationOperator})
    if length(a) == 1
        comm(a[1], b)
    elseif length(a) == 2
        CompositeTerm(a[1]) * comm(a[2], b) +
        comm(a[1], b) * CompositeTerm(a[2])
    else
        sum(
            prod(CompositeTerm, a[1:i-1]; init=CompositeTerm(1)) *
            comm(a[i], b) *
            prod(CompositeTerm, a[i+1:end]; init=CompositeTerm(1))
            for i in eachindex(a)
        )
    end
end

function comm(a::CompositeTerm{A}, b::CompositeTerm{B}) where
{A<:Number,B<:Number}
    if isempty(a.operators) || isempty(b.operators)
        zero(promote_type(A, B))
    else
        a_inds = get_all_inds(a)
        b_inds = get_all_inds(b)

        a_sum_inds_ovlp = intersect(a.sum_inds, b_inds)
        b_sum_inds_ovlp = intersect(b.sum_inds, a_inds)

        ex_table_a = [i => ind(i.o, i.n * "₁") for i in a_sum_inds_ovlp]
        ex_table_b = [i => ind(i.o, i.n * "₂") for i in b_sum_inds_ovlp]

        a = exchange_index(a, ex_table_a)
        b = exchange_index(b, ex_table_b)

        nonop = get_nonop(a) * get_nonop(b)
        op = comm(a.operators, b.operators)
        summation(get_nosum(nonop) * op, nonop.sum_inds)
    end
end

function comm(a::CompositeTerm{A}, b::SumType{B}) where {A<:Number,B<:Number}
    sum(comm(a, t) for t in b.terms)
end

function comm(a::SumType{A}, b::CompositeTerm{B}) where {A<:Number,B<:Number}
    sum(comm(t, b) for t in a.terms)
end

function comm(a::SumType{A}, b::SumType{B}) where {A<:Number,B<:Number}
    sum(comm(t1, t2) for t1 in a.terms, t2 in b.terms)
end
