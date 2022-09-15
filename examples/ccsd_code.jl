using PyCall
using LinearAlgebra: norm
using OMEinsum
using LinearAlgebra: I
np = pyimport("numpy")
pyscf = pyimport("pyscf")

struct OEI_int
    oo
    ov
    vo
    vv
end

struct ERI_int
    oooo
    ooov
    oovo
    oovv
    ovoo
    ovov
    ovvo
    ovvv
    vooo
    voov
    vovo
    vovv
    vvoo
    vvov
    vvvo
    vvvv
end

# ENERGY
function energy(t_vovo, F, g, h)
    u_vovo = 2 * t_vovo - permutedims(t_vovo, (1,4,3,2));
    E_ = [mol_h2o.energy_nuc()];
    E_ += np.einsum("iajb,aibj->", g.ovov, u_vovo, optimize="optimal") * 1.0;
    E_ += np.einsum("ii->", F.oo, optimize="optimal") * 1.0;
    E_ += np.einsum("ii->", h.oo, optimize="optimal") * 1.0;
    return E_;
end

# OMEGA SINGLES
function singles(t_vovo, F, g)
    u_vovo = 2 * t_vovo - permutedims(t_vovo, (1,4,3,2));
    omega_ai = np.einsum("ai->ai", F.vo, optimize="optimal") * 1.0;
    omega_ai += np.einsum("abjc,bicj->ai", g.vvov, u_vovo, optimize="optimal") * 1.0;
    omega_ai += np.einsum("jb,aibj->ai", F.ov, u_vovo, optimize="optimal") * 1.0;
    omega_ai += np.einsum("jbki,akbj->ai", g.ovoo, u_vovo, optimize="optimal") * -1.0;
    return norm(omega_ai);
end

# OMEGA DOUBLES
function doubles(t_vovo, F, g)
    omega_aibj = np.einsum("aibj->aibj", g.vovo, optimize="optimal") * 1.0;
    omega_aibj += np.einsum("ac,bjci->aibj", F.vv, t_vovo, optimize="optimal") * 1.0;
    omega_aibj += np.einsum("bc,aicj->aibj", F.vv, t_vovo, optimize="optimal") * 1.0;
    omega_aibj += np.einsum("acbd,cidj->aibj", g.vvvv, t_vovo, optimize="optimal") * 1.0;
    omega_aibj += np.einsum("kcld,aibk,cjdl->aibj", g.ovov, t_vovo, t_vovo, optimize="optimal") * -2.0;
    omega_aibj += np.einsum("kcld,aibk,cldj->aibj", g.ovov, t_vovo, t_vovo, optimize="optimal") * 1.0;
    omega_aibj += np.einsum("kcld,aicj,bkdl->aibj", g.ovov, t_vovo, t_vovo, optimize="optimal") * -2.0;
    omega_aibj += np.einsum("kcld,aicj,bldk->aibj", g.ovov, t_vovo, t_vovo, optimize="optimal") * 1.0;
    omega_aibj += np.einsum("kcld,aick,bjdl->aibj", g.ovov, t_vovo, t_vovo, optimize="optimal") * 4.0;
    omega_aibj += np.einsum("kcld,aick,bldj->aibj", g.ovov, t_vovo, t_vovo, optimize="optimal") * -2.0;
    omega_aibj += np.einsum("kcld,aicl,bjdk->aibj", g.ovov, t_vovo, t_vovo, optimize="optimal") * -2.0;
    omega_aibj += np.einsum("kcld,aicl,bkdj->aibj", g.ovov, t_vovo, t_vovo, optimize="optimal") * 1.0;
    omega_aibj += np.einsum("kcld,akbj,cidl->aibj", g.ovov, t_vovo, t_vovo, optimize="optimal") * -2.0;
    omega_aibj += np.einsum("kcld,akbj,cldi->aibj", g.ovov, t_vovo, t_vovo, optimize="optimal") * 1.0;
    omega_aibj += np.einsum("kcld,akbl,cidj->aibj", g.ovov, t_vovo, t_vovo, optimize="optimal") * 1.0;
    omega_aibj += np.einsum("kcld,akci,bjdl->aibj", g.ovov, t_vovo, t_vovo, optimize="optimal") * -2.0;
    omega_aibj += np.einsum("kcld,akci,bldj->aibj", g.ovov, t_vovo, t_vovo, optimize="optimal") * 1.0;
    omega_aibj += np.einsum("kcld,akcl,bjdi->aibj", g.ovov, t_vovo, t_vovo, optimize="optimal") * 1.0;
    omega_aibj += np.einsum("kcld,akdi,bjcl->aibj", g.ovov, t_vovo, t_vovo, optimize="optimal") * 1.0;
    omega_aibj += np.einsum("kcld,akdj,blci->aibj", g.ovov, t_vovo, t_vovo, optimize="optimal") * 1.0;
    omega_aibj += np.einsum("kcld,akdl,bjci->aibj", g.ovov, t_vovo, t_vovo, optimize="optimal") * -2.0;
    omega_aibj += np.einsum("acki,bjck->aibj", g.vvoo, t_vovo, optimize="optimal") * -1.0;
    omega_aibj += np.einsum("ackj,bkci->aibj", g.vvoo, t_vovo, optimize="optimal") * -1.0;
    omega_aibj += np.einsum("aikc,bjck->aibj", g.voov, t_vovo, optimize="optimal") * 2.0;
    omega_aibj += np.einsum("aikc,bkcj->aibj", g.voov, t_vovo, optimize="optimal") * -1.0;
    omega_aibj += np.einsum("bcki,akcj->aibj", g.vvoo, t_vovo, optimize="optimal") * -1.0;
    omega_aibj += np.einsum("bckj,aick->aibj", g.vvoo, t_vovo, optimize="optimal") * -1.0;
    omega_aibj += np.einsum("bjkc,aick->aibj", g.voov, t_vovo, optimize="optimal") * 2.0;
    omega_aibj += np.einsum("bjkc,akci->aibj", g.voov, t_vovo, optimize="optimal") * -1.0;
    omega_aibj += np.einsum("ki,akbj->aibj", F.oo, t_vovo, optimize="optimal") * -1.0;
    omega_aibj += np.einsum("kj,aibk->aibj", F.oo, t_vovo, optimize="optimal") * -1.0;
    omega_aibj += np.einsum("kilj,akbl->aibj", g.oooo, t_vovo, optimize="optimal") * 1.0;
    for a = 1:nv, i=1:no
        omega_aibj[a,i,a,i] /= 2
    end
    return norm(omega_aibj);
end

function t1_transform(h, t1)
    nv = size(t1, 1)
    no = size(t1, 2)
    nmo = no + nv
    t = zeros(nmo, nmo)
    t[no+1:end, 1:no] .= t1
    x = I - t
    y = I + t'
    ht1 = x * h * y'
    return ht1
end

function t1_transform_eri(g, t1)
    nv = size(t1, 1)
    no = size(t1, 2)
    nmo = no + nv
    t = zeros(nmo, nmo)
    t[no+1:end, 1:no] .= t1
    x = I - t
    y = I + t'
    @ein gt1[p,q,r,s] := (((g[t,u,m,n] * x[p,t]) * y[q,u]) * x[r,m]) * y[s,n]
    return gt1
end

function t1_fock(t1, h, g)
    no = size(t1, 2)
    Ft1 = t1_transform(h, t1)
    gt1 = t1_transform_eri(g, t1)
    for p = axes(Ft1,1), q = axes(Ft1,2)
        Ft1[p,q] += sum(2*gt1[p,q,i,i] - gt1[p,i,i,q] for i = 1:no)
    end
    return Ft1
end

println("PYSCF CCSD")
mol_h2o = pyscf.gto.M(atom = "O 0 0 0; H 0 1 0; H 0 0 1", basis = "ccpvdz");
rhf_h2o = pyscf.scf.RHF(mol_h2o);
e_h2o = rhf_h2o.kernel();
mycc = rhf_h2o.CCSD().run()

t_vo = permutedims(mycc.t1, (2,1))
t_vovo = permutedims(mycc.t2, (3,1,4,2));

eri = mol_h2o.intor("int2e");
eris = pyscf.ao2mo.kernel(eri, rhf_h2o.mo_coeff);
eri_t1 = t1_transform_eri(eris, t_vo);

h_ao = mol_h2o.intor_symmetric("int1e_kin") + mol_h2o.intor_symmetric("int1e_nuc");
hmo = np.einsum("pi,pq,qj->ij", rhf_h2o.mo_coeff, h_ao, rhf_h2o.mo_coeff);
ht1 = t1_transform(hmo, t_vo);
Ft1 = t1_fock(t_vo, hmo, eris)

no = mol_h2o.nelectron ÷ 2
nv = size(hmo, 1) - no

ht1_oei = OEI_int(
    ht1[1:no,1:no],
    ht1[1:no,no+1:end],
    ht1[no+1:end,1:no],
    ht1[no+1:end,no+1:end]
);

Ft1_oei = OEI_int(
    Ft1[1:no,1:no],
    Ft1[1:no,no+1:end],
    Ft1[no+1:end,1:no],
    Ft1[no+1:end,no+1:end]
);

gt1 = ERI_int(
    eri_t1[    1:no,     1:no,     1:no,     1:no],
    eri_t1[    1:no,     1:no,     1:no, no+1:end],
    eri_t1[    1:no,     1:no, no+1:end,     1:no],
    eri_t1[    1:no,     1:no, no+1:end, no+1:end],
    eri_t1[    1:no, no+1:end,     1:no,     1:no],
    eri_t1[    1:no, no+1:end,     1:no, no+1:end],
    eri_t1[    1:no, no+1:end, no+1:end,     1:no],
    eri_t1[    1:no, no+1:end, no+1:end, no+1:end],
    eri_t1[no+1:end,     1:no,     1:no,     1:no],
    eri_t1[no+1:end,     1:no,     1:no, no+1:end],
    eri_t1[no+1:end,     1:no, no+1:end,     1:no],
    eri_t1[no+1:end,     1:no, no+1:end, no+1:end],
    eri_t1[no+1:end, no+1:end,     1:no,     1:no],
    eri_t1[no+1:end, no+1:end,     1:no, no+1:end],
    eri_t1[no+1:end, no+1:end, no+1:end,     1:no],
    eri_t1[no+1:end, no+1:end, no+1:end, no+1:end]
);

println("\nExcitationOperators CCSD")
@show energy(t_vovo, Ft1_oei, gt1, ht1_oei);
@show energy(t_vovo, Ft1_oei, gt1, ht1_oei) - [e_h2o];
@show singles(t_vovo, Ft1_oei, gt1);
@show doubles(t_vovo, Ft1_oei, gt1);