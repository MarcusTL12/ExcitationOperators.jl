using Test

using ExcitationOperators

println("INDICES:")

i = ind(occ, "i")
j = ind(occ, "j")
a = ind(vir, "a")
p = ind(gen, "p")

@show i, j, a, p
@show i < a
@show i > a
@show i <= a
@show i >= a
@show i > j

println("\nEXCITATION OPERATORS:")

Eai = E(a, i)
Epj = E(p, j)

@show Eai, Epj

println("\nKROENECKER DELTA:")

d1 = δ(p, i)
d2 = δ(i, p)
d3 = δ(a, i)

@show d1, d2, d3
@show d1 == d2
@show d1 > d2
@show d2 >= d1

println("\nPRODUCTS:")

dd = d1 * d2

@show dd

dd *= Epj

@show dd

dd *= Eai

@show dd

@show dd'

@show 3dd

dd2 = d2 * E(i, j)

@show dd2

@show E(i, j) * 0

println("\nSUMS:")

dd3 = dd - dd2

@show dd3

dd4 = dd3 + E(a, p)

@show dd4

dd5 = δ(p, a) - E(j, i)

dd6 = nothing + dd5 + 3dd4 + nothing - 3E(a, p)

@show dd6
@show -dd6' + 2//3

println("\nCOMMUTATORS:")

q = ind(gen, "q")
b = ind(vir, "b")

@show comm(E(p, q), E(a, i))
@show comm(E(i, a), E(b, j))
@show comm(E(i, a), E(j, b))
@show comm(E(i, a), E(a, j))
