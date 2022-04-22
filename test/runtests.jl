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
Epj = E(p, j)

@show Eai, Epj

d1 = δ(p, i)
d2 = δ(i, p)
d3 = δ(a, i)

@show d1, d2, d3
@show d1 == d2
@show d1 > d2
@show d2 <= d3

dd = d1 * d3

@show dd

dd *= Epj

@show dd

dd *= Eai

@show dd

@show dd'

@show 3dd

dd2 = d2 * E(i, j)

@show dd2

dd3 = dd - dd2

@show dd3

dd4 = dd3 + E(a, p)

@show dd4

dd5 = δ(p, a) - E(j, i)

dd6 = nothing + dd5 + 3dd4 + nothing - 3E(a, p)

@show dd6
@show -dd6' + 2//3
