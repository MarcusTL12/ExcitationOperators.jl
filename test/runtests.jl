using Test

using ExcitationOperators

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

Eai = E(a, i)

@show Eai

d1 = δ(p, i)
d2 = δ(i, p)
d3 = δ(a, i)

@show d1, d2, d3
@show d1 == d2
@show d1 > d2
@show d2 <= d3

dd = d1 * d3

@show dd
