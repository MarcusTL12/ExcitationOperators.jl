export printlatex, printlnlatex, latexstring

function printlatex(io::IO, i::MOIndex, color)
    if color && !isgen(i)
        print(io, "{\\color{", isocc(i) ? "green" : "cyan", '}', i.n, '}')

    else
        print(io, i.n)
    end
end

function printlatex(io::IO, t::Tensor, color)
    print(io, get_symbol(t), "_{")
    inds = get_indices(t)
    printlatex(io, first(inds), color)
    for i in inds[2:end]
        print(io, ' ')
        printlatex(io, i, color)
    end
    print(io, '}')
end

function printlatex(io::IO, d::KroeneckerDelta, color)
    print(io, "\\delta_{")
    printlatex(io, d.p, color)
    printlatex(io, d.q, color)
    print(io, '}')
end

function printlatex(io::IO, o::ExcitationOperator, color)
    print(io, "E_{")
    printlatex(io, o.p, color)
    printlatex(io, o.q, color)
    print(io, '}')
end

function printlatex(io::IO, s::T) where {T<:Number}
    print(io, s)
end

function printlatex(io::IO, s::Rational{T}) where {T}
    if isone(denominator(s))
        print(io, numerator(s))
    else
        print(io, "\\frac{", numerator(s), "}{", denominator(s), '}')
    end
end

function printlatex(io::IO, t::CompositeTerm{T}, color=false) where {T<:Number}
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
            printlatex(io, t.scalar)
        end
    elseif all_nonscalar_empty
        print(io, t.scalar)
    end

    if !isempty(t.sum_inds)
        printsep()
        print(io, "\\sum_{")
        inds = collect(t.sum_inds)
        printlatex(io, first(inds), color)
        for i in inds[2:end]
            print(io, ' ')
            printlatex(io, i, color)
        end
        print(io, "}{")
        sep[] = false
    end

    for d in t.deltas
        printsep()
        printlatex(io, d, color)
    end

    for ten in t.tensors
        printsep()
        printlatex(io, ten, color)
    end

    for op in t.operators
        printsep()
        printlatex(io, op, color)
    end

    if !isempty(t.sum_inds)
        print(io, '}')
    end
end

function printlatex(io::IO, s::SumType{T}, color=false) where {T<:Number}
    printlatex(io, first(s.terms), color)
    for t in s.terms[2:end]
        if t.scalar < zero(T)
            print(io, "\n- ")
            printlatex(io, -t, color)
        else
            print(io, "\n+ ")
            printlatex(io, t, color)
        end
    end
end

function printlatex(v, color=false)
    printlatex(stdout, v, color)
end

function printlnlatex(io::IO, v, color=false)
    printlatex(io, v, color)
    print(io, '\n')
end

function printlnlatex(v, color=false)
    printlnlatex(stdout, v, color)
end

function latexstring(v, color=false)
    buf = IOBuffer()
    printlatex(buf, v, color)
    String(take!(buf))
end
