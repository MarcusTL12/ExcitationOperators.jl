struct CompositeTerm{T<:Number}
    scalar::T
    sum_inds::SortedSet{MOIndex}
    deltas::SortedSet{KroeneckerDelta}
    tensors::Vector{Tensor}
    operators::Vector{ExcitationOperator}

    function CompositeTerm(
        n::T,
        sum_inds::SortedSet{MOIndex},
        deltas::SortedSet{KroeneckerDelta},
        tensors::Vector{Tensor},
        operators::Vector{ExcitationOperator}) where {T<:Number}
        if iszero(n)
            new{T}(
                n,
                SortedSet{MOIndex}(),
                SortedSet{KroeneckerDelta}(),
                Tensor[],
                ExcitationOperator[]
            )
        else
            delta_inds = Set{MOIndex}()
            for d in deltas
                push!(delta_inds, d.p)
                push!(delta_inds, d.q)
            end

            t = new{T}(n, sum_inds, deltas, sort(tensors), operators)

            to_sum_over = nothing
            for i in sum_inds
                if i ∈ delta_inds
                    to_sum_over = i
                    break
                end
            end

            if isnothing(to_sum_over)
                t
            else
                delete!(t.sum_inds, to_sum_over)
                t = summation(t, to_sum_over)

                t
            end
        end
    end
end

function CompositeTerm(
    n::T,
    deltas::SortedSet{KroeneckerDelta},
    tensors::Vector{Tensor},
    operators::Vector{ExcitationOperator}) where {T<:Number}
    CompositeTerm(
        n,
        SortedSet{MOIndex}(),
        deltas,
        tensors,
        operators
    )
end

function Base.:(==)(a::CompositeTerm{A}, b::CompositeTerm{B}) where
{A<:Number,B<:Number}
    (a.sum_inds, a.scalar, a.deltas, a.tensors, a.operators) ==
    (b.sum_inds, b.scalar, b.deltas, b.tensors, b.operators)
end

function Base.isless(a::CompositeTerm{A}, b::CompositeTerm{B}) where
{A<:Number,B<:Number}
    (a.sum_inds, a.deltas, a.tensors, a.operators, -a.scalar) <
    (b.sum_inds, b.deltas, b.tensors, b.operators, -b.scalar)
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

    if !isempty(t.sum_inds)
        printsep()
        print(io, "∑_")
        for i in t.sum_inds
            print(io, i.n)
        end
        print(io, '(')
        sep[] = false
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

    if !isempty(t.sum_inds)
        print(io, ')')
    end
end

Base.iszero(t::CompositeTerm{T}) where {T<:Number} = iszero(t.scalar)

Base.adjoint(t::CompositeTerm{T}) where {T<:Number} = CompositeTerm(
    t.scalar',
    t.sum_inds,
    t.deltas,
    Tensor[adjoint(ten) for ten in t.tensors],
    reverse(adjoint.(t.operators))
)

# Chech whether the non-scalar part is the same
function issimilar(a::CompositeTerm{A}, b::CompositeTerm{B}) where
{A<:Number,B<:Number}
    (a.sum_inds, a.deltas, a.tensors, a.operators) ==
    (b.sum_inds, b.deltas, b.tensors, b.operators)
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

CompositeTerm(operators::Vector{ExcitationOperator}) = CompositeTerm(
    1,
    SortedSet{KroeneckerDelta}(),
    Tensor[],
    operators
)

# Make constructors for types to make CompositeTerm the external interface type

export E, δ, real_tensor

E(p::MOIndex, q::MOIndex) = CompositeTerm(ExcitationOperator(p, q))
δ(p::MOIndex, q::MOIndex) = CompositeTerm(KroeneckerDelta(p, q))
real_tensor(symbol, indices...) =
    CompositeTerm(RealTensor(symbol, collect(indices)))

# Overloading multiplication

function Base.:*(a::A, b::CompositeTerm{B}) where {A<:Number,B<:Number}
    CompositeTerm(a * b.scalar, b.sum_inds, b.deltas, b.tensors, b.operators)
end

function Base.:*(a::CompositeTerm{A}, b::B) where {A<:Number,B<:Number}
    b * a
end

function Base.:*(a::CompositeTerm{A}, b::CompositeTerm{B}) where
{A<:Number,B<:Number}
    CompositeTerm(
        a.scalar * b.scalar,
        union(a.sum_inds, b.sum_inds),
        union(a.deltas, b.deltas),
        Tensor[a.tensors; b.tensors],
        [a.operators; b.operators]
    )
end

# Negating scalar

Base.:-(t::CompositeTerm{T}) where {T<:Number} =
    CompositeTerm(-t.scalar, t.sum_inds, t.deltas, t.tensors, t.operators)

# Utility method for promoting scalar type of term. Will just try calling
# NT(scalar)
function convert_scalar(::Type{NT}, t::CompositeTerm{T}) where
{NT<:Number,T<:Number}
    CompositeTerm(NT(t.scalar), t.sum_inds, t.deltas, t.tensors, t.operators)
end

# Utility method for getting the non-operator part of a term
function get_nonop(t::CompositeTerm{T}) where {T<:Number}
    CompositeTerm(
        t.scalar,
        t.sum_inds,
        t.deltas,
        t.tensors,
        ExcitationOperator[])
end

# Exchanging indices

function exchange_index(t::CompositeTerm{T}, from::MOIndex, to::MOIndex) where
{T<:Number}
    sum_inds = SortedSet([i == from ? to : i for i in t.sum_inds])

    deltas = SortedSet{KroeneckerDelta}()
    for d in t.deltas
        nd = exchange_index(d, from, to)
        if nd isa KroeneckerDelta
            push!(deltas, nd)
        elseif iszero(nd)
            return CompositeTerm(zero(A))
        end
    end

    tensors = Tensor[exchange_index(ten, from, to) for ten in t.tensors]
    operators = [exchange_index(o, from, to) for o in t.operators]

    CompositeTerm(t.scalar, sum_inds, deltas, tensors, operators)
end

function exchange_index(t::CompositeTerm{T}, mapping) where {T<:Number}
    foldl((acc, (from, to)) -> exchange_index(acc, from, to), mapping; init=t)
end
