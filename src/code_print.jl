export print_code

function print_code(t :: CompositeTerm, symbol :: String, to :: String)
    scalar_str = float(t.scalar)

    # Make einsum_str
    einsum_str = "\""
    for a in t.tensors
        for b in get_indices(a)
            einsum_str *= string(b.n)
        end
        einsum_str *= ","
    end
    einsum_str = einsum_str[begin:end-1] * "->" * to * "\""

    # Make tensor_str
    tensor_str = ""
    for a in t.tensors
        tensor_str *= ", $(get_symbol(a))_"
        for b in get_indices(a)
            tensor_str *= string(b.n)
        end
    end

    tensor_str = replace(tensor_str, "a" => "v", "b" => "v", "c" => "v", "d" => "v", "i" => "o", "j" => "o", "k" => "o", "l" => "o")
    tensor_str = replace(tensor_str, "F_" => "F.", "h_" => "h.", "g_" => "g.")
    println("$(symbol)_$(to) += np.einsum($einsum_str$tensor_str, optimize=\"optimal\") * $scalar_str;")
end
