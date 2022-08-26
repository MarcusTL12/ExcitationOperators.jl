export comm

function comm(::A, ::B) where {A<:Number,B}
    zero(promote_type(A, B))
end

function comm(::A, ::B) where {A,B<:Number}
    zero(promote_type(A, B))
end

function comm(a::ExcitationOperator, b::ExcitationOperator)
    mul_collide(δ(a.q, b.p), E(a.p, b.q)) -
    mul_collide(δ(a.p, b.q), E(b.p, a.q))
end

function comm(a::ExcitationOperator, b::Vector{ExcitationOperator})
    if length(b) == 1
        comm(a, b[1])
    elseif length(b) == 2
        mul_collide(CompositeTerm(b[1]), comm(a, b[2])) +
        mul_collide(comm(a, b[1]), CompositeTerm(b[2]))
    else
        acc = CompositeTerm(0)

        for i in eachindex(b)
            left = CompositeTerm(1)
            for t in b[1:i-1]
                left = mul_collide(left, CompositeTerm(t))
            end

            right = CompositeTerm(1)
            for t in b[1:i-1]
                right = mul_collide(right, CompositeTerm(t))
            end

            mid = comm(a, b[i])

            acc += mul_collide(left, mul_collide(mid, right))
        end

        acc
    end
end

function comm(a::Vector{ExcitationOperator}, b::Vector{ExcitationOperator})
    if length(a) == 1
        comm(a[1], b)
    elseif length(a) == 2
        mul_collide(CompositeTerm(a[1]), comm(a[2], b)) +
        mul_collide(comm(a[1], b), CompositeTerm(a[2]))
    else
        acc = CompositeTerm(0)

        for i in eachindex(a)
            left = CompositeTerm(1)
            for t in a[1:i-1]
                left = mul_collide(left, CompositeTerm(t))
            end

            right = CompositeTerm(1)
            for t in a[1:i-1]
                right = mul_collide(right, CompositeTerm(t))
            end

            mid = comm(a[i], b)

            acc += mul_collide(left, mul_collide(mid, right))
        end

        acc
    end
end

function comm(a::CompositeTerm{A}, b::CompositeTerm{B}) where
{A<:Number,B<:Number}
    if isempty(a.operators) || isempty(b.operators)
        zero(promote_type(A, B))
    else
        nonop = mul_collide(get_nonop(a), get_nonop(b))
        op = comm(a.operators, b.operators)
        mul_collide(nonop, op)
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
