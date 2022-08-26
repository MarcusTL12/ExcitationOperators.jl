using Test
using ExcitationOperators
using ExcitationOperators.BasicStuff.StandardIndices
using ExcitationOperators.BasicStuff.StandardOperators

ei = ind(vir, "e")
f = ind(vir, "f")
m = ind(occ, "m")
n = ind(occ, "n")
o1 = ind(occ, "o₁")
v1 = ind(vir, "v₁")
o2 = ind(occ, "o₂")
v2 = ind(vir, "v₂")

tbar = real_tensor("tbar", a, i)
tbar2 = sym_tensor("tbar", a, i, b, j)
T2  = summation(1//2 * sym_tensor("t", c, k, d, l) * E(c, k) * E(d, l), [c, k, d, l])
T22 = summation(1//2 * sym_tensor("t", ei, m, f, n) * E(ei, m) * E(f, n), [ei, m, f, n])

P1 = 1//2 * E(i, a)
P2 = (1//3 * E(i, a) * E(j, b) + 1//6 * E(i, b) * E(j, a)) 

Λ = 1 + summation(P1 * tbar, [a,i]) + 1//2 * summation(P2 * tbar2, [a,i,b,j])
BCH2(A, B, C) = A + comm(A, B) + 1//2 * comm(comm(A, B), C)

for pp in [o1, v1], qq in [o2, v2]
    Epq = Λ * BCH2(E(pp,qq), T2, T22)
    D = exval(Epq)
    D = D |> simplify |> sort_sum_sym_tensor
    println("D_$(pp.n)$(qq.n) = $D")
end