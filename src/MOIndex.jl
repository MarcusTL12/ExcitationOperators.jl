export occ, vir, gen, ind

@enum Occupation vir gen occ

struct MOIndex
    o::Occupation
    n::String
end

ind(o, name) = MOIndex(o, name)

Base.isless(p::MOIndex, q::MOIndex) = (p.o, p.n) < (q.o, q.n)
