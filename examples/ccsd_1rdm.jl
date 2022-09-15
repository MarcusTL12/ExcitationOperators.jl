using Test
using ExcitationOperators
using ExcitationOperators.BasicStuff.StandardIndices
using ExcitationOperators.BasicStuff.StandardOperators

o1 = ind(occ, "o₁")
v1 = ind(vir, "v₁")
o2 = ind(occ, "o₂")
v2 = ind(vir, "v₂")

a1 = ind(vir, "a₁")
i1 = ind(occ, "i₁")
b1 = ind(vir, "b₁")
j1 = ind(occ, "j₁")

a2 = ind(vir, "a₂")
i2 = ind(occ, "i₂")
b2 = ind(vir, "b₂")
j2 = ind(occ, "j₂")

tbar = cc_amp_tensor("tbar", a, i)
tbar2 = cc_amp_tensor("tbar", a, i, b, j)
T2 = summation(1//2 * cc_amp_tensor("t", a1, i1, b1, j1) * E(a1, i1) * E(b1, j1), [a1, i1, b1, j1])
T22 = summation(1//2 * cc_amp_tensor("t", a2, i2, b2, j2) * E(a2, i2) * E(b2, j2), [a2, i2, b2, j2])

P1 = 1//2 * E(i, a)
P2 = (1//3 * E(i, a) * E(j, b) + 1//6 * E(i, b) * E(j, a)) 

Λ = 1 + summation(P1 * tbar, [a,i]) + 1//2 * summation(P2 * tbar2, [a,i,b,j])
BCH2(A, B, C) = A + comm(A, B) + 1//2 * comm(comm(A, B), C)

for pp in [o1, v1], qq in [o2, v2]
    Epq = Λ * BCH2(E(pp,qq), T2, T22)
    D = exval(Epq)
    D = D |> simplify
    println("D_$(pp.n)$(qq.n) = $D")
end