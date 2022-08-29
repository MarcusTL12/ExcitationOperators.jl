export occ, vir, gen, ind, make_occ, make_vir

@enum Occupation gen vir occ

struct MOIndex
    o::Occupation
    n::String
end

ind(o, name) = MOIndex(o, name)

Base.isless(p::MOIndex, q::MOIndex) = (p.o, p.n) < (q.o, q.n)

#

function Base.show(io::IO, i::MOIndex)
    if isvir(i)
        print(io, "\x1b[36m")
    elseif isocc(i)
        print(io, "\x1b[92m")
    end
    print(io, i.n, "\x1b[39m")
end

# Utility method for promoting general index to occupied or virtual index

function make_occ(i::MOIndex)
    if i.o == gen
        ind(occ, i.n)
    else
        i
    end
end

function make_vir(i::MOIndex)
    if i.o == gen
        ind(vir, i.n)
    else
        i
    end
end

# TODO: refactor code to use these
isgen(i::MOIndex) = i.o == gen
isvir(i::MOIndex) = i.o == vir
isocc(i::MOIndex) = i.o == occ
