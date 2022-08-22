struct SumType{T<:Number}
    terms::Vector{CompositeTerm{T}}

    function SumType(terms::Vector{CompositeTerm{T}}) where {T<:Number}
        counter = Dict{
            Tuple{
                SortedSet{KroeneckerDelta},
                Vector{Tensor},
                Vector{ExcitationOperator}
            },
            T
        }()

        for t in terms
            k = (t.deltas, t.tensors, t.operators)
            counter[k] = get(counter, k, zero(T)) + t.scalar
        end

        terms = [
            CompositeTerm(s, d, t, o)
            for ((d, t, o), s) in counter if !iszero(s)
        ]
        sort!(terms)

        if isempty(terms)
            0
        elseif length(terms) == 1
            terms[1]
        else
            new{T}(terms)
        end
    end
end

function Base.show(io::IO, s::SumType{T}) where {T<:Number}
    print(io, first(s.terms))
    for t in s.terms[2:end]
        if t.scalar < zero(T)
            print(io, " - ", -t)
        else
            print(io, " + ", t)
        end
    end
end

# Overloading +

function Base.:+(a::CompositeTerm{A}, b::CompositeTerm{B}) where
{A<:Number,B<:Number}
    T = promote_type(A, B)
    SumType([convert_scalar(T, a), convert_scalar(T, b)])
end

function Base.:+(a::SumType{A}, b::CompositeTerm{B}) where {A<:Number,B<:Number}
    T = promote_type(A, B)
    SumType([convert_scalar.(T, a.terms); convert_scalar(T, b)])
end

function Base.:+(a::CompositeTerm{A}, b::SumType{B}) where {A<:Number,B<:Number}
    b + a
end

function Base.:+(a::SumType{A}, b::SumType{B}) where {A<:Number,B<:Number}
    T = promote_type(A, B)
    SumType([convert_scalar.(T, a.terms); convert_scalar.(T, b.terms)])
end

function Base.:+(a::Union{CompositeTerm{A},SumType{A}}, b::B) where
{A<:Number,B<:Number}
    a + CompositeTerm(b)
end

function Base.:+(a::A, b::Union{CompositeTerm{B},SumType{B}}) where
{A<:Number,B<:Number}
    b + a
end

# Overloading -

function Base.:-(a::SumType{A}) where {A<:Number}
    SumType([-t for t in a.terms])
end

function Base.:-(
    a::Union{A,CompositeTerm{A},SumType{A}},
    b::Union{B,CompositeTerm{B},SumType{B}}) where {A<:Number,B<:Number}
    a + (-b)
end

# Overloading *

function Base.:*(a::Union{A,CompositeTerm{A}}, b::SumType{B}) where
{A<:Number,B<:Number}
    SumType([a * t for t in b.terms])
end

function Base.:*(a::SumType{A}, b::Union{B,CompositeTerm{B}}) where
{A<:Number,B<:Number}
    b * a
end

function Base.:*(a::SumType{A}, b::SumType{B}) where {A<:Number,B<:Number}
    sum(t1 * t2 for t1 in a.terms, t2 in b.terms)
end

# Adjoint

function Base.adjoint(a::SumType{A}) where {A<:Number}
    SumType(adjoint.(a.terms))
end
