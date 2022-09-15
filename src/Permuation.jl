export expand_perm, exchange_index_safe, collapse_perm

# Safe as long as no indices with names similar to iᵗ already exist
function exchange_index_safe(t, mapping)
    tmp_indices = [ind(q.o, q.n * "ᵗ") for (_, q) in mapping]
    tmp_mapping = [p => q for ((p, _), q) in zip(mapping, tmp_indices)]
    final_mapping = [p => q for (p, (_, q)) in zip(tmp_indices, mapping)]

    exchange_index(exchange_index(t, tmp_mapping), final_mapping)
end

function expand_perm(
    t::CompositeTerm{T},
    perm_inds::Vector{Tuple{MOIndex,MOIndex}}
) where {T<:Number}
    sum(begin
        mapping = Pair{MOIndex,MOIndex}[]
        for i in eachindex(perm_inds)
            if perm(i) != i
                push!(mapping, perm_inds[i][1] => perm_inds[perm(i)][1])
                push!(mapping, perm_inds[i][2] => perm_inds[perm(i)][2])
            end
        end
        exchange_index_safe(t, mapping)
    end for perm in PermGen(length(perm_inds)))
end

function collapse_perm_first_term(
    s::SumType{T},
    perm_inds::Vector{Tuple{MOIndex,MOIndex}}
) where {T<:Number}
    perm_inds_line = MOIndex[]
    for (p, q) in perm_inds
        push!(perm_inds_line, p)
        push!(perm_inds_line, q)
    end
    permop = cc_amp_tensor("P", perm_inds_line...)

    lookup = Set(s.terms)

    first_term = first(s.terms)
    expanded_term = simplify(expand_perm(first_term, perm_inds))
    if expanded_term isa SumType &&
       length(expanded_term.terms) == factorial(length(perm_inds))
        if all(t ∈ lookup for t in expanded_term.terms)
            remaining_terms = setdiff(s.terms, expanded_term.terms)
            return (true, permop * first_term, SumType(remaining_terms))
        end
    end

    (false, first_term, SumType(s.terms[2:end]))
end

function collapse_perm(
    s::SumType{T},
    perm_inds::Vector{Tuple{MOIndex,MOIndex}}
) where {T<:Number}
    rest = s
    acc = CompositeTerm(0)

    while rest isa SumType
        _, contracted, new_rest = collapse_perm_first_term(rest, perm_inds)
        rest = new_rest
        acc += contracted
    end

    acc + rest
end
