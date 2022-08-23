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
        t = if iszero(n)
            new{T}(
                n,
                SortedSet{MOIndex}(),
                SortedSet{KroeneckerDelta}(),
                Tensor[],
                ExcitationOperator[]
            )
        else
            t = new{T}(n, sum_inds, deltas, sort(tensors), operators)

            delta_inds = Set{MOIndex}()
            for d in t.deltas
                push!(delta_inds, d.p)
                push!(delta_inds, d.q)
            end

            to_sum_over = nothing
            for i in t.sum_inds
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
    (collect(a.sum_inds), collect(a.deltas), a.tensors, a.operators, -a.scalar) <
    (collect(b.sum_inds), collect(b.deltas), b.tensors, b.operators, -b.scalar)
end

function Base.show(io::IO, t::CompositeTerm{T}) where {T<:Number}
    sep = Ref(false)

    function printsep()
        if sep[]
            print(io, ' ')
        end
        sep[] = true
    end

    all_nonscalar_empty = isempty(t.sum_inds) && isempty(t.deltas) &&
                          isempty(t.tensors) && isempty(t.operators)

    if !isone(t.scalar)
        if isone(-t.scalar)
            print(io, '-')
        else
            printsep()
            print(io, t.scalar)
        end
    elseif all_nonscalar_empty
        print(io, t.scalar)
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

Base.zero(::Type{CompositeTerm{T}}) where {T<:Number} = CompositeTerm(zero(T))

# Make constructors for types to make CompositeTerm the external interface type

export E, δ, real_tensor, e

E(p::MOIndex, q::MOIndex) = CompositeTerm(ExcitationOperator(p, q))
δ(p::MOIndex, q::MOIndex) = CompositeTerm(KroeneckerDelta(p, q))
real_tensor(symbol, indices...) =
    CompositeTerm(RealTensor(symbol, collect(indices)))

e(p::MOIndex, q::MOIndex, r::MOIndex, s::MOIndex) =
    E(p, q) * E(r, s) - δ(q, r) * E(p, s)

# Overloading multiplication

function Base.:*(a::A, b::CompositeTerm{B}) where {A<:Number,B<:Number}
    CompositeTerm(a * b.scalar, b.sum_inds, b.deltas, b.tensors, b.operators)
end

function Base.:*(a::CompositeTerm{A}, b::B) where {A<:Number,B<:Number}
    b * a
end

function Base.:*(a::CompositeTerm{A}, b::CompositeTerm{B}) where
{A<:Number,B<:Number}
    common_sum_inds = intersect(a.sum_inds, b.sum_inds)

    ex_table_a = [i => ind(i.o, i.n * "₁") for i in common_sum_inds]
    ex_table_b = [i => ind(i.o, i.n * "₂") for i in common_sum_inds]

    a = exchange_index(a, ex_table_a)
    b = exchange_index(b, ex_table_b)

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

# Index cleanup

function check_general_indices(t::CompositeTerm{T}) where {T<:Number}
    exchange_table = Pair{MOIndex,MOIndex}[]
    for d in t.deltas
        if d.p.o == gen && d.q.o != gen
            new_ind = if d.q.o == occ
                make_occ(d.p)
            else
                make_vir(d.p)
            end
            push!(exchange_table, d.p => new_ind)
        elseif d.p.o != gen && d.q.o == gen
            new_ind = if d.p.o == occ
                make_occ(d.q)
            else
                make_vir(d.q)
            end
            push!(exchange_table, d.q => new_ind)
        end
    end

    exchange_index(t, exchange_table)
end

export cleanup_indices

function cleanup_indices(
    t::CompositeTerm{T};
    gen_queue=["p", "q", "r", "s", "t", "u", "v", "w"],
    occ_queue=["i", "j", "k", "l", "m", "n"],
    vir_queue=["a", "b", "c", "d", "e", "f"]
) where {T<:Number}
    t = check_general_indices(t)

    gen_queue = [ind(gen, n) for n in gen_queue]
    occ_queue = [ind(occ, n) for n in occ_queue]
    vir_queue = [ind(vir, n) for n in vir_queue]

    nonsum_inds = Set{MOIndex}()

    for d in t.deltas
        push!(nonsum_inds, d.p)
        push!(nonsum_inds, d.q)
    end

    for t in t.tensors
        for i in get_indices(t)
            push!(nonsum_inds, i)
        end
    end

    for o in t.operators
        push!(nonsum_inds, o.p)
        push!(nonsum_inds, o.q)
    end

    setdiff!(nonsum_inds, t.sum_inds)

    setdiff!(gen_queue, nonsum_inds)
    setdiff!(occ_queue, nonsum_inds)
    setdiff!(vir_queue, nonsum_inds)

    tmp_ex_table = [i => ind(i.o, i.n * "t") for i in t.sum_inds]

    t = exchange_index(t, tmp_ex_table)

    ex_table = Pair{MOIndex,MOIndex}[]

    for (_, i) in tmp_ex_table
        new_ind = if i.o == gen
            popfirst!(gen_queue)
        elseif i.o == occ
            popfirst!(occ_queue)
        else
            popfirst!(vir_queue)
        end

        push!(ex_table, i => new_ind)
    end

    exchange_index(t, ex_table)
end

# Experimental splitting for simplification. Splits sum over general indices
# into one sum over occupied and one over virtual.
export split_summation

function split_summation(t::CompositeTerm{T}) where {T<:Number}
    gen_ind = nothing
    for i in t.sum_inds
        if i.o == gen
            gen_ind = i
            break
        end
    end

    if isnothing(gen_ind)
        t
    else
        t_occ = exchange_index(t, gen_ind, make_occ(gen_ind))
        t_vir = exchange_index(t, gen_ind, make_vir(gen_ind))

        split_summation(t_occ) + split_summation(t_vir)
    end
end


# Experimental combining of sums to simplify. Tries to identify whether the two
# terms are similar except for one summation index that is occupied in one sum
# and virtual in the other, then combines them to a single sum
export combine_summation

subscript(i) = join(Char(0x2080 + d) for d in reverse!(digits(i)))

function make_comb_ind_func()
    counter = 0
    function make_comb_ind()
        counter += 1
        ind(gen, "ζ" * subscript(counter))
    end
end

const make_comb_ind = make_comb_ind_func()

function combine_summation(a::CompositeTerm{A}, b::CompositeTerm{B}) where
{A<:Number,B<:Number}
    if a.scalar == b.scalar
        a_has_occ = any(isocc(i) for i in a.sum_inds)
        a_has_vir = any(isvir(i) for i in a.sum_inds)
        b_has_occ = any(isocc(i) for i in b.sum_inds)
        b_has_vir = any(isvir(i) for i in b.sum_inds)

        test_ind = ind(gen, "tₑₛₜ")

        if a_has_occ && b_has_vir || a_has_vir && b_has_occ
            for i in a.sum_inds, j in b.sum_inds
                if isocc(i) && isvir(j) || isvir(i) && isocc(j)
                    a_test = cleanup_indices(exchange_index(a, i, test_ind))
                    b_test = cleanup_indices(exchange_index(b, j, test_ind))

                    if a_test == b_test
                        return cleanup_indices(
                            exchange_index(a, i, make_comb_ind())
                        )
                    end
                end
            end
        end
    end
end

function combine_summation(a::CompositeTerm{A}) where {A<:Number}
    a
end
