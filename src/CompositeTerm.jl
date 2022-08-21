struct CompositeTerm{T<:Number}
    scalar::T
    deltas::SortedSet{KroeneckerDelta}
    tensors::Vector{Tensor}
    operators::Vector{ExcitationOperator}

    function CompositeTerm(n::T, deltas::SortedSet{KroeneckerDelta},
        tensors::Vector{Tensor}, operators::Vector{ExcitationOperator}) where
    {T<:Number}
        if iszero(n)
            new{T}(
                n,
                SortedSet{KroeneckerDelta}(),
                Tensor[],
                ExcitationOperator[]
            )
        else
            new{T}(n, deltas, sort(tensors), operators)
        end
    end
end

function Base.:(==)(a::CompositeTerm{A}, b::CompositeTerm{B}) where
{A<:Number,B<:Number}
    (a.scalar, a.deltas, a.tensors, a.operators) ==
    (b.scalar, b.deltas, b.tensors, b.operators)
end

function Base.isless(a::CompositeTerm{A}, b::CompositeTerm{B}) where
{A<:Number,B<:Number}
    (a.scalar, a.deltas, a.tensors, a.operators) <
    (b.scalar, b.deltas, b.tensors, b.operators)
end

function Base.show(io::IO, t::CompositeTerm{T}) where {T<:Number}
    sep = Ref(false)

    function printsep()
        if sep[]
            print(io, ' ')
        end
        sep[] = true
    end

    if !isone(t.scalar)
        if isone(-t.scalar)
            print(io, '-')
        else
            printsep()
            print(io, t.scalar)
        end
    end

    for d in t.deltas
        printsep()
        print(io, d)
    end

    for ten in t.tensors
        printsep()
        print(io, ten)
    end

    for op in t.operators
        printsep()
        print(io, op)
    end
end

Base.iszero(t::CompositeTerm{T}) where {T<:Number} = iszero(t.scalar)

Base.adjoint(t::CompositeTerm{T}) where {T<:Number} = CompositeTerm(
    t.scalar',
    t.deltas,
    Tensor[adjoint(ten) for ten in t.tensors],
    reverse(adjoint.(t.operators))
)

# Chech whether the non-scalar part is the same
function issimilar(a::CompositeTerm{A}, b::CompositeTerm{B}) where
{A<:Number,B<:Number}
    (a.deltas, a.tensors, a.operators) ==
    (b.deltas, b.tensors, b.operators)
end

# Promotation

CompositeTerm(n::T) where {T<:Number} = CompositeTerm(
    n,
    SortedSet{KroeneckerDelta}(),
    Tensor[],
    ExcitationOperator[]
)

CompositeTerm(d::KroeneckerDelta) = CompositeTerm(
    1,
    SortedSet{KroeneckerDelta}([d]),
    Tensor[],
    ExcitationOperator[]
)

CompositeTerm(t::T) where {T<:Tensor} = CompositeTerm(
    1,
    SortedSet{KroeneckerDelta}(),
    Tensor[t],
    ExcitationOperator[]
)

CompositeTerm(e::ExcitationOperator) = CompositeTerm(
    1,
    SortedSet{KroeneckerDelta}(),
    Tensor[],
    [e]
)

# Make constructors for types to make CompositeTerm the external interface type

export E, δ, real_tensor

E(p::MOIndex, q::MOIndex) = CompositeTerm(ExcitationOperator(p, q))
δ(p::MOIndex, q::MOIndex) = CompositeTerm(KroeneckerDelta(p, q))
real_tensor(symbol, indices...) =
    CompositeTerm(RealTensor(symbol, collect(indices)))

# Overloading multiplication

function Base.:*(a::A, b::CompositeTerm{B}) where {A<:Number,B<:Number}
    CompositeTerm(a * b.scalar, b.deltas, b.tensors, b.operators)
end

function Base.:*(a::CompositeTerm{A}, b::B) where {A<:Number,B<:Number}
    b * a
end

function Base.:*(a::CompositeTerm{A}, b::CompositeTerm{B}) where
{A<:Number,B<:Number}
    CompositeTerm(
        a.scalar * b.scalar,
        union(a.deltas, b.deltas),
        Tensor[a.tensors; b.tensors],
        [a.operators; b.operators]
    )
end

# Negating scalar

Base.:-(t::CompositeTerm{T}) where {T<:Number} =
    CompositeTerm(-t.scalar, t.deltas, t.tensors, t.operators)
