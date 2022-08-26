using Test
using ExcitationOperators
using ExcitationOperators.BasicStuff.StandardIndices
using ExcitationOperators.BasicStuff.StandardOperators

using ExcitationOperators: CompositeTerm, SumType, MOIndex, exchange_index
function aibj_symmetry(sumtype :: Union{CompositeTerm{A},SumType{A}}, a :: MOIndex, i :: MOIndex, b :: MOIndex, j :: MOIndex) where {A<:Number}
    # Implements the symmetry ai <-> bj
    bt = ind(vir, "bt")
    jt = ind(vir, "jt")
    count = 1
    while count ≤ length(sumtype.terms)
        len = length(sumtype.terms)
        sumtype.terms[count] = ExcitationOperators.exchange_index(sumtype.terms[count], a, bt)
        sumtype.terms[count] = ExcitationOperators.exchange_index(sumtype.terms[count], b, a)
        sumtype.terms[count] = ExcitationOperators.exchange_index(sumtype.terms[count], bt, b)
        sumtype.terms[count] = ExcitationOperators.exchange_index(sumtype.terms[count], i, jt)
        sumtype.terms[count] = ExcitationOperators.exchange_index(sumtype.terms[count], j, i)
        sumtype.terms[count] = ExcitationOperators.exchange_index(sumtype.terms[count], jt, j)
        sumtype = sum(sumtype.terms) |> simplify |> sort_sum_sym_tensor
        if typeof(sumtype) <: SumType
            if length(sumtype.terms) < len
                len = length(sumtype.terms)
            else
                count += 1
            end
        else
            break
        end
    end
    return sumtype
end

ei = ind(vir, "e")
f = ind(vir, "f")
m = ind(occ, "m")
n = ind(occ, "n")
o = ind(occ, "o")

taibj = sym_tensor("t", c, k, d, l)
T2 = summation(1//2 * taibj * E(c, k) * E(d, l), [c, k, d, l])
taibj2 = sym_tensor("t", ei, m, f, n)
T22 = summation(1//2 * taibj2 * E(ei, m) * E(f, n), [ei, m, f, n])

hF = summation( (real_tensor("F", p, q) + summation(-2//1 * sym_tensor("g", p,q,o,o) + 1//1 * sym_tensor("g", p,o,o,q), [o])) * E(p, q), [p,q])
HF = hF + g

# Correct
T2u = summation(1//2 * (2//3 * sym_tensor("u", c, k, d, l) + 1//3 * sym_tensor("u", c, l, d, k)) * E(c, k) * E(d, l), [c, k, d, l])
ECCSD = exval( 1//2 * H + 1//2 * HF + comm(HF, T2u)) |> simplify |> sort_sum_sym_tensor
@show ECCSD;

# Correct
T2u = summation(1//2 * (2//3 * sym_tensor("u", c, k, d, l) + 1//3 * sym_tensor("u", c, l, d, k)) * E(c, k) * E(d, l), [c, k, d, l])
P1 = 1//2 * E(i, a)
Ω1 = exval( P1 * (HF + comm(HF, T2u))) |> simplify |> sort_sum_sym_tensor
@show Ω1;

# Correct
P2 = (1//3 * E(i, a) * E(j, b) + 1//6 * E(i, b) * E(j, a)) 
Ω2 = exval( P2 * (HF + comm(HF, T2) + 1//2 * comm(comm(HF, T2), T22))) |> simplify |> sort_sum_sym_tensor
Ω2 = aibj_symmetry(Ω2, a, i, b, j)
@show Ω2;


# Pinkbook (Helgaker CCSD expressions)
# Verify that ExcitationOperators.jl produce reasonable expressions for CCSD
ECCSD_pink = ( summation(real_tensor("h", i, i) + real_tensor("F", i, i), i) + summation(sym_tensor("g", i, a, j, b) * sym_tensor("u", a, i, b, j), [a, i, b, j]) ) |> simplify |> sort_sum_sym_tensor
@test (ECCSD == ECCSD_pink);

ΩA1 = summation(sym_tensor("u", c, k, d, i) * sym_tensor("g", a, d, k, c), [c, k, d])
ΩB1 = summation(-sym_tensor("u", a, k, c, l) * sym_tensor("g", k, i, l, c), [c, k, l])
ΩC1 = summation(sym_tensor("u", a, i, c, k) * real_tensor("F", k, c), [k, c])
ΩD1 = real_tensor("F", a, i)
Ω1_pink = (ΩA1 + ΩB1 + ΩC1 + ΩD1) |> simplify |> sort_sum_sym_tensor
@test (Ω1 == Ω1_pink);

ΩA2 = sym_tensor("g", a, i, b, j) + summation( sym_tensor("t", c, i, d, j) * sym_tensor("g", a, c, b, d), [c, d])
ΩB2 = summation( sym_tensor("t", a, k, b, l) * ( sym_tensor("g", k, i, l, j) + summation(sym_tensor("t", c, i, d, j) * sym_tensor("g", k, c, l, d), [c, d])), [k, l])
ΩC2 = - 1//2 * summation(sym_tensor("t", b, k, c, j) * (sym_tensor("g", k, i, a, c) - 1//2 * summation(sym_tensor("t", a, l, d, i) * sym_tensor("g", k, d, l, c), [d, l])), [c, k]) -
               summation(sym_tensor("t", b, k, c, i) * (sym_tensor("g", k, j, a, c) - 1//2 * summation(sym_tensor("t", a, l, d, j) * sym_tensor("g", k, d, l, c), [d, l])), [c, k])
ubjck = 2//1 * sym_tensor("t", b, j, c, k) - sym_tensor("t", b, k, c, j)
uaidl = 2//1 * sym_tensor("t", a, i, d, l) - sym_tensor("t", a, l, d, i)
Laikc = 2//1 * sym_tensor("g", a, i, k, c) - sym_tensor("g", a, c, k, i)
Lldkc = 2//1 * sym_tensor("g", l, d, k, c) - sym_tensor("g", l, c, k, d)
ΩD2 = 1//2 * summation(ubjck * (Laikc + 1//2 * summation(uaidl * Lldkc, [d, l])), [c, k])

ubkdl = 2//1 * sym_tensor("t", b, k, d, l) - sym_tensor("t", b, l, d, k)
ucldj = 2//1 * sym_tensor("t", c, l, d, j) - sym_tensor("t", c, j, d, l)
ΩE2 = summation(sym_tensor("t", a, i, c, j) * (real_tensor("F", b, c) - summation(ubkdl * sym_tensor("g", l, d, k, c), [d, k, l])), [c]) - 
      summation(sym_tensor("t", a, i, b, k) * (real_tensor("F", k, j) + summation(ucldj * sym_tensor("g", k, d, l, c), [c, d, l])), [k])

function Paibj(sumterms, a, i, b, j)
    # Exchanges ai <-> bj
    sumterms2 = deepcopy(sumterms)
    bt = ind(vir, "bt")
    jt = ind(vir, "jt")
    for index = 1:length(sumterms2.terms)
        sumterms2.terms[index] = ExcitationOperators.exchange_index(sumterms2.terms[index], a, bt)
        sumterms2.terms[index] = ExcitationOperators.exchange_index(sumterms2.terms[index], b, a)
        sumterms2.terms[index] = ExcitationOperators.exchange_index(sumterms2.terms[index], bt, b)
        sumterms2.terms[index] = ExcitationOperators.exchange_index(sumterms2.terms[index], i, jt)
        sumterms2.terms[index] = ExcitationOperators.exchange_index(sumterms2.terms[index], j, i)
        sumterms2.terms[index] = ExcitationOperators.exchange_index(sumterms2.terms[index], jt, j)
    end
    return sumterms + sumterms2
end

Ω2_pink = ΩA2 + ΩB2 + Paibj(ΩC2 + ΩD2 + ΩE2, a, i, b, j)
Ω2_pink = aibj_symmetry(Ω2_pink, a, i, b, j) |> simplify |> sort_sum_sym_tensor
@test (Ω2 == Ω2_pink);