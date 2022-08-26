export cleanup_indices, split_summation, combine_summation, simplify,
    sort_sum_sym_tensor

function cleanup_indices(
    t::CompositeTerm{T};
    gen_queue=["p", "q", "r", "s", "t", "u", "v", "w"],
    occ_queue=["i", "j", "k", "l", "m", "n", "o", "π"],
    vir_queue=["a", "b", "c", "d", "e", "f", "g", "h"]
) where {T<:Number}
    t = check_general_indices(t)

    gen_queue = [ind(gen, n) for n in gen_queue]
    occ_queue = [ind(occ, n) for n in occ_queue]
    vir_queue = [ind(vir, n) for n in vir_queue]

    nonsum_inds = Set{MOIndex}()
    sorted_sum_inds = MOIndex[]

    function add_sorted_sum_ind(i)
        if i ∈ t.sum_inds && i ∉ sorted_sum_inds
            push!(sorted_sum_inds, i)
        end
    end

    for d in t.deltas
        push!(nonsum_inds, d.p)
        push!(nonsum_inds, d.q)
        add_sorted_sum_ind(d.p)
        add_sorted_sum_ind(d.q)
    end

    for t in t.tensors
        for i in get_indices(t)
            push!(nonsum_inds, i)
            add_sorted_sum_ind(i)
        end
    end

    for o in t.operators
        push!(nonsum_inds, o.p)
        push!(nonsum_inds, o.q)
        add_sorted_sum_ind(o.p)
        add_sorted_sum_ind(o.q)
    end

    for i in t.sum_inds
        add_sorted_sum_ind(i)
    end

    setdiff!(nonsum_inds, t.sum_inds)

    setdiff!(gen_queue, nonsum_inds)
    setdiff!(occ_queue, nonsum_inds)
    setdiff!(vir_queue, nonsum_inds)

    tmp_ex_table = [i => ind(i.o, i.n * "t") for i in sorted_sum_inds]

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

function cleanup_indices(
    s::SumType{T};
    gen_queue=["p", "q", "r", "s", "t", "u", "v", "w"],
    occ_queue=["i", "j", "k", "l", "m", "n"],
    vir_queue=["a", "b", "c", "d", "e", "f"]
) where {T<:Number}
    sum(cleanup_indices(
        t;
        gen_queue=copy(gen_queue),
        occ_queue=copy(occ_queue),
        vir_queue=copy(vir_queue)
    ) for t in s.terms)
end

# Experimental splitting for simplification. Splits sum over general indices
# into one sum over occupied and one over virtual.
function split_summation(s::SumType{T}) where {T<:Number}
    sum(split_summation(t) for t in s.terms)
end

function combine_summation_single_pass(a::SumType{T}) where {T<:Number}
    terms = copy(a.terms)

    for i in eachindex(terms), j in i+1:length(terms)
        comb = combine_summation(terms[i], terms[j])
        if !isnothing(comb)
            deleteat!(terms, j)
            deleteat!(terms, i)
            return (true, comb + sum(terms))
        end
    end

    (false, a)
end

combine_summation_single_pass(a::CompositeTerm{T}) where {T<:Number} =
    (false, a)

function combine_summation(a::SumType{T}) where {T<:Number}
    done, comb = combine_summation_single_pass(a)
    while done
        done, comb = combine_summation_single_pass(comb)
    end
    comb
end

function simplify(a::Union{CompositeTerm{A},SumType{A}}) where {A<:Number}
    a |> cleanup_indices |> split_summation |>
    cleanup_indices |> combine_summation
end

# Attempt at simplifying
# ∑_aij g_iajj = ∑_aij g_iija
function sort_sum_sym_tensor(t::CompositeTerm{T}) where {T<:Number}
    sum_inds_gen = MOIndex[i for i in t.sum_inds if isgen(i)]
    sum_inds_occ = MOIndex[i for i in t.sum_inds if isocc(i)]
    sum_inds_vir = MOIndex[i for i in t.sum_inds if isvir(i)]

    tmp_inds_gen = [ind(i.o, i.n * "ₜ") for i in sum_inds_gen]
    tmp_inds_occ = [ind(i.o, i.n * "ₜ") for i in sum_inds_occ]
    tmp_inds_vir = [ind(i.o, i.n * "ₜ") for i in sum_inds_vir]

    tmp_ex_table = [i => ind(i.o, i.n * "ₜ") for i in t.sum_inds]

    make_ex_table(tmp_inds, sum_inds, perm) = [
        tmp_inds[i] => sum_inds[perm(i)]
        for i in eachindex(sum_inds)
    ]

    t_tmp = exchange_index(t, tmp_ex_table)

    for gen_perm in PermGen(length(sum_inds_gen)),
        occ_perm in PermGen(length(sum_inds_occ)),
        vir_perm in PermGen(length(sum_inds_vir))

        ex_table = [
            make_ex_table(tmp_inds_gen, sum_inds_gen, gen_perm)
            make_ex_table(tmp_inds_occ, sum_inds_occ, occ_perm)
            make_ex_table(tmp_inds_vir, sum_inds_vir, vir_perm)
        ]

        t_perm = exchange_index(t_tmp, ex_table)

        if t_perm < t
            t = t_perm
        end
    end

    t
end

function sort_sum_sym_tensor(s::SumType{T}) where {T<:Number}
    sum(sort_sum_sym_tensor(t) for t in s.terms)
end

