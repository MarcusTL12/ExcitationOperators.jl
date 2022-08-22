export occ, vir, gen, ind

@enum Occupation gen vir occ

struct MOIndex
    o::Occupation
    n::String
end

ind(o, name) = MOIndex(o, name)

Base.isless(p::MOIndex, q::MOIndex) = (p.o, p.n) < (q.o, q.n)

# Utility method for promoting general index to occupied or virtual index

function make_occ(i::MOIndex)
    if i.o == gen
        ind(occ, i.n * "ᵒ")
    else
        i
    end
end

function make_vir(i::MOIndex)
    if i.o == gen
        ind(vir, i.n * "ᵛ")
    else
        i
    end
end

# TODO: make is_occ, is_gen and is_occ methods
# and refactor code to make use of them
