using ExcitationOperators
using ExcitationOperators.BasicStuff.StandardIndices
using ExcitationOperators.BasicStuff.StandardOperators

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

# Complicated
# Missing ai <-> bj symmetry in Ω2
P2 = (1//3 * E(i, a) * E(j, b) + 1//6 * E(i, b) * E(j, a)) 
exval(P2 * HF) |> simplify
exval(P2 * comm(HF, T2)) |> simplify |> sort_sum_sym_tensor
exval(P2 * 1//2 * comm(comm(HF, T2), T22)) |> simplify |> sort_sum_sym_tensor
Ω2 = exval( P2 * (HF + comm(HF, T2) + 1//2 * comm(comm(HF, T2), T22))) |> simplify |> sort_sum_sym_tensor
@show Ω2;