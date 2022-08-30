module BasicStuff

module StandardIndices
using ExcitationOperators

export p, q, r, s, i, j, k, l, a, b, c, d

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
end

module StandardOperators
using ExcitationOperators
using ..StandardIndices

export ht, gt, h, g, H

# Tensors
ht(p, q) = sym_tensor("h", p, q)
gt(p, q, r, s) = sym_tensor("g", p, q, r, s)

# Operators
h = summation(ht(p, q) * E(p, q), [p, q])
g = summation(1 // 2 * gt(p, q, r, s) * e(p, q, r, s), [p, q, r, s])

H = h + g

end

end