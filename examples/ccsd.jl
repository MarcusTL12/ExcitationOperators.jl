using Test
using ExcitationOperators
using ExcitationOperators.BasicStuff.StandardIndices
using ExcitationOperators.BasicStuff.StandardOperators

a1 = ind(vir, "a₁")
i1 = ind(occ, "i₁")
b1 = ind(vir, "b₁")
j1 = ind(occ, "j₁")

a2 = ind(vir, "a₂")
i2 = ind(occ, "i₂")
b2 = ind(vir, "b₂")
j2 = ind(occ, "j₂")

o = ind(occ, "o")

T2 = summation(1//2 * cc_amp_tensor("t", a1, i1, b1, j1) * E(a1, i1) * E(b1, j1), [a1, i1, b1, j1])
T22 = summation(1//2 * cc_amp_tensor("t", a2, i2, b2, j2) * E(a2, i2) * E(b2, j2), [a2, i2, b2, j2])
T2u = summation(1//2 * (2//3 * cc_amp_tensor("u", a1, i1, b1, j1) + 1//3 * cc_amp_tensor("u", a1, j1, b1, i1)) * E(a1, i1) * E(b1, j1), [a1, i1, b1, j1])

hF = summation( (real_tensor("F", p, q) + summation(-2//1 * psym_tensor("g", p,q,o,o) + psym_tensor("g", p,o,o,q), [o])) * E(p, q), [p,q])
HF = hF + g

# Correct
ECCSD = exval( 1//2 * H + 1//2 * HF + comm(HF, T2u)) |> simplify
@show ECCSD;

# Correct
P1 = 1//2 * E(i, a)
Ω1 = exval( P1 * (HF + comm(HF, T2u))) |> simplify
@show Ω1;

# Correct
P2 = (1//3 * E(i, a) * E(j, b) + 1//6 * E(i, b) * E(j, a)) 
Ω2 = exval( P2 * (HF + comm(HF, T2) + 1//2 * comm(comm(HF, T2), T22))) |> simplify
@show Ω2;

# Pinkbook (Helgaker CCSD expressions)
# Verify that ExcitationOperators.jl produce reasonable expressions for CCSD
function Paibj(sumterms)
    # P_aibj A_abij = A_abij + A_baji
    sumterms2 = deepcopy(sumterms)
    for index = 1:length(sumterms2.terms)
        sumterms2.terms[index] = exchange_index_safe(sumterms2.terms[index], [a=>b, b=>a, i=>j, j=>i])
    end
    return sumterms + sumterms2
end

ECCSD_pink = ( summation(real_tensor("h", i, i) + real_tensor("F", i, i), i) + summation(psym_tensor("g", i, a, j, b) * psym_tensor("u", a, i, b, j), [a, i, b, j]) ) |> simplify
@test (ECCSD == ECCSD_pink);

ΩA1 = summation(cc_amp_tensor("u", c, k, d, i) * psym_tensor("g", a, d, k, c), [c, k, d])
ΩB1 = summation(-cc_amp_tensor("u", a, k, c, l) * psym_tensor("g", k, i, l, c), [c, k, l])
ΩC1 = summation(cc_amp_tensor("u", a, i, c, k) * real_tensor("F", k, c), [k, c])
ΩD1 = real_tensor("F", a, i)
Ω1_pink = (ΩA1 + ΩB1 + ΩC1 + ΩD1) |> simplify
@test (Ω1 == Ω1_pink);

ΩA2 = psym_tensor("g", a, i, b, j) + summation( cc_amp_tensor("t", c, i, d, j) * psym_tensor("g", a, c, b, d), [c, d])
ΩB2 = summation( cc_amp_tensor("t", a, k, b, l) * ( psym_tensor("g", k, i, l, j) + summation(cc_amp_tensor("t", c, i, d, j) * psym_tensor("g", k, c, l, d), [c, d])), [k, l])
ΩC2 = - 1//2 * summation(cc_amp_tensor("t", b, k, c, j) * (psym_tensor("g", k, i, a, c) - 1//2 * summation(cc_amp_tensor("t", a, l, d, i) * psym_tensor("g", k, d, l, c), [d, l])), [c, k]) -
               summation(cc_amp_tensor("t", b, k, c, i) * (psym_tensor("g", k, j, a, c) - 1//2 * summation(psym_tensor("t", a, l, d, j) * psym_tensor("g", k, d, l, c), [d, l])), [c, k])
ubjck = 2//1 * psym_tensor("t", b, j, c, k) - psym_tensor("t", b, k, c, j)
uaidl = 2//1 * psym_tensor("t", a, i, d, l) - psym_tensor("t", a, l, d, i)
Laikc = 2//1 * psym_tensor("g", a, i, k, c) - psym_tensor("g", a, c, k, i)
Lldkc = 2//1 * psym_tensor("g", l, d, k, c) - psym_tensor("g", l, c, k, d)
ΩD2 = 1//2 * summation(ubjck * (Laikc + 1//2 * summation(uaidl * Lldkc, [d, l])), [c, k])
ubkdl = 2//1 * psym_tensor("t", b, k, d, l) - psym_tensor("t", b, l, d, k)
ucldj = 2//1 * psym_tensor("t", c, l, d, j) - psym_tensor("t", c, j, d, l)
ΩE2 = summation(psym_tensor("t", a, i, c, j) * (real_tensor("F", b, c) - summation(ubkdl * psym_tensor("g", l, d, k, c), [d, k, l])), [c]) - 
      summation(psym_tensor("t", a, i, b, k) * (real_tensor("F", k, j) + summation(ucldj * psym_tensor("g", k, d, l, c), [c, d, l])), [k])
Ω2_pink = ΩA2 + ΩB2 + Paibj(ΩC2 + ΩD2 + ΩE2) |> simplify
@test (Ω2 == Ω2_pink);

# Introduce Paibj back into the equations
PΩ2 = collapse_perm(Ω2, [(a,i),(b,j)])
@show PΩ2;

# Code generation
for t in ECCSD.terms
    print_code(t, "E", "")
end

for t in Ω1.terms
    print_code(t, "omega", "ai")
end

for t in Ω2.terms
    print_code(t, "omega", "aibj")
end