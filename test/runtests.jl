using Test

using ExcitationOperators

println("INDICES:")

@show Int.([vir, gen, occ])

p = ind(gen, "p")
q = ind(gen, "q")
r = ind(gen, "r")
s = ind(gen, "s")
i = ind(occ, "i")
j = ind(occ, "j")
k = ind(occ, "k")
l = ind(occ, "l")
a = ind(vir, "a")
b = ind(vir, "b")
c = ind(vir, "c")
d = ind(vir, "d")

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

println("\nTENSORS:")

hpq = real_tensor("h", p, q)
hrs = real_tensor("h", r, s)
gpqrs = real_tensor("g", p, q, r, s)

@show hpq hrs gpqrs

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

@show E(p, q) * 3

@show E(p, q) * hpq
@show E(p, q) * E(r, s) * gpqrs
@show -E(p, q) * E(r, s) * δ(q, r) * gpqrs
@show hpq * gpqrs
@show hrs * gpqrs

println("\nSUMS:")

@show E(i, j) + 0

dd3 = dd - dd2

@show dd3

dd4 = dd3 + E(a, p)

@show dd4

dd5 = gpqrs * δ(p, a) - E(j, i) * hpq

@show dd5
@show -dd5' + 2 // 3

@show (3 + δ(i, j) * E(a, b) - 2 // 3 * E(p, q)) * (3 * δ(p, i) - E(q, j))

println("\nCOMMUTATORS:")

@show comm(E(p, q), E(a, i))
@show comm(E(i, a), E(b, j))
@show comm(E(i, a), E(j, b))
@show comm(E(i, a), E(a, j))

println()

@show comm(E(i, a), E(p, q) * E(b, j))
@show comm(E(i, a) * E(p, q), E(b, j))
@show comm(E(i, a) * E(p, q), E(a, i) * E(b, j))

println()

@show comm(E(i, a), comm(E(p, q), E(b, j)))
@show comm(comm(E(i, a), E(p, q)), E(b, j))

println()

@show comm(E(i, a) * E(j, b), comm(E(p, q), E(c, k) * E(d, l)))

println("\nSIMPLE SUMMATIONS:")

@show summation(hpq * δ(p, q) * E(p, q), p)

println()
@show tmp = hpq * comm(E(i, a), E(p, q) * E(b, j))
@show summation(tmp, Set([p, q]))

# println("\nEXPECTATION VALUES:")
# @show exval(E(p, q))
# @show exval(E(a, q))
# @show exval(E(p, i))

# println()

# ex1 = exval(δ(a, b) * E(p, q) * E(r, s))

# @show ex1

# ex2 = ex1 * exval(E(p, q))

# @show ex2

# @show exval(E(p, q)) * exval(E(p, q))
# @show exval(E(r, s)) * exval(E(p, q))

# @show exval(E(p, q) * Tensor2("h", p, q))
# @show exval(E(p, q) * E(r, s) * Tensor4("g", p, q, r, s))

# println("\nEXPECTATION SUMS:")

# @show -ex1
# @show 3ex2

# @show ex1 - 3
# @show ex2 - ex1

# @show ex2 - δ(i, j)

# println()

# @show exval(comm(E(i, a), comm(E(p, q), E(b, j))))

# println()

# @show exval(comm(E(i, a) * E(j, b), comm(E(p, q), E(c, k) * E(d, l))))