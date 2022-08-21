struct SumType{T<:Number}
    terms::Vector{CompositeTerm{T}}

    function SumType(terms::Vector{CompositeTerm{T}}) where {T<:Number}
        new{T}(sort(terms))
    end
end

# Overloading + and -

