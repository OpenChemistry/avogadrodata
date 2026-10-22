import os, sys, itertools
import numpy as np
from rdkit import Chem
from rdkit.Geometry import Point3D
from spiral import *

out = sys.argv[1]; os.makedirs(out, exist_ok=True)

def write(name, symbols, X, bonds, orders=None, title=None):
    m = Chem.RWMol()
    for s in symbols: m.AddAtom(Chem.Atom(s))
    for k, (a, b) in enumerate(bonds):
        o = 1 if orders is None else orders[k]
        m.AddBond(int(a), int(b), {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE}[o])
    conf = Chem.Conformer(len(symbols))
    for i, p in enumerate(X): conf.SetAtomPosition(i, Point3D(*[float(c) for c in p]))
    m.AddConformer(conf)
    m.SetProp("_Name", name)
    open(os.path.join(out, name + ".mol"), "w").write(Chem.MolToV3KMolBlock(m.GetMol()))
    print(name, len(symbols), "atoms", len(bonds), "bonds", flush=True)

def faces_of_atoms(faces, nvd):
    """for each dual vertex, the list of atoms (triangles) around it"""
    around = [[] for _ in range(nvd)]
    for t, f in enumerate(faces):
        for v in f: around[v].append(t)
    return around

def pentagon_bonds(bonds, around, degs):
    pb = set()
    for v, d in enumerate(degs):
        if d == 5:
            A = set(around[v])
            for (a, b) in bonds:
                if a in A and b in A: pb.add((a, b))
    return pb

def build(n, ipr_search):
    (adj, faces), pent = ipr_search
    nv, bonds, faces2 = fullerene_graph((adj, faces))
    assert nv == n and len(bonds) == 3*n//2
    degs = [len(a) for a in adj]
    around = faces_of_atoms(faces, len(adj))
    return adj, faces, bonds, degs, around

# ---- C60, C70 ----
res = {}
for n, name in ((60, "c60"), (70, "c70")):
    s = search_ipr2(n)
    adj, faces, bonds, degs, around = build(n, s)
    pb = pentagon_bonds(bonds, around, degs)
    target = {b: (1.45 if b in pb else 1.39) for b in bonds}
    X = embed(n, bonds, target, seed=1)
    assert X is not None, name
    res[n] = (adj, faces, bonds, degs, around, X)
    orders = [1 if b in pb else 2 for b in bonds] if n == 60 else None
    if n == 60: assert sum(orders) - len(orders) == 30
    write(name, ["C"]*n, X, bonds, orders)
    L = [np.linalg.norm(X[a]-X[b]) for a, b in bonds]
    print("  bond lengths", min(L), max(L))

# ---- C240, C540 by subdividing the C60 dual triangulation ----
adj, faces, bonds60, degs, around, X60 = res[60]
nd = len(adj)
P = np.zeros((nd, 3))
for v in range(nd):
    P[v] = X60[around[v]].mean(axis=0)
P /= np.linalg.norm(P, axis=1)[:, None]

def subdivide(k):
    pos = {}
    tris = []
    def key(tid, w):
        nz = tuple(sorted((f, c) for f, c in zip(faces[tid], w) if c))
        if len(nz) == 3: nz = nz + (tid,)
        return nz
    def xyz(tid, w):
        return sum(c*P[f] for f, c in zip(faces[tid], w))/k
    for tid in range(len(faces)):
        def pt(a, b):
            w = (k-a-b, a, b)
            kk = key(tid, w)
            if kk not in pos: pos[kk] = xyz(tid, w)
            return kk
        for a in range(k):
            for b in range(k-a):
                tris.append((pt(a, b), pt(a+1, b), pt(a, b+1)))
                if a+b <= k-2:
                    tris.append((pt(a+1, b), pt(a, b+1), pt(a+1, b+1)))
    edge = {}
    for t, tr in enumerate(tris):
        for a, b in itertools.combinations(tr, 2):
            edge.setdefault(frozenset((a, b)), []).append(t)
    bonds = sorted({tuple(sorted(v)) for v in edge.values() if len(v) == 2})
    cent = np.array([sum(pos[v] for v in tr) for tr in tris])
    cent /= np.linalg.norm(cent, axis=1)[:, None]
    return bonds, cent

for k, name in ((2, "c240"), (3, "c540")):
    bonds, cent = subdivide(k)
    nvert = len(cent)
    R = 1.42*np.sqrt(nvert/(4*np.pi)*2.6)
    # relax from the spherical start (single try: start is already unfolded)
    nb = [[] for _ in range(nvert)]
    for a, b in bonds: nb[a].append(b); nb[b].append(a)
    assert set(len(x) for x in nb) == {3}, name
    target = {b: 1.42 for b in bonds}
    # use embed's energy through a one-try call with an explicit start
    import spiral
    X = None
    rng = np.random.default_rng(0)
    # monkeypatch: run embed with custom initial points
    def embed_from(X0):
        bset = set(bonds)
        pairs13 = set()
        for v in range(nvert):
            for a, b in itertools.combinations(nb[v], 2): pairs13.add((min(a, b), max(a, b)))
        b_idx = np.array(bonds); p13 = np.array(sorted(pairs13))
        allpairs = np.array([(i, j) for i in range(nvert) for j in range(i+1, nvert)
                             if (i, j) not in bset and (i, j) not in pairs13])
        def energy(x):
            Xc = x.reshape(-1, 3); g = np.zeros_like(Xc); E = 0.0
            d = Xc[b_idx[:, 0]]-Xc[b_idx[:, 1]]; r = np.linalg.norm(d, axis=1); dd = r-1.42
            E += 10*(dd**2).sum(); gg = (20*dd/r)[:, None]*d
            np.add.at(g, b_idx[:, 0], gg); np.add.at(g, b_idx[:, 1], -gg)
            d = Xc[p13[:, 0]]-Xc[p13[:, 1]]; r = np.linalg.norm(d, axis=1); dd = r-2.40
            E += 2*(dd**2).sum(); gg = (4*dd/r)[:, None]*d
            np.add.at(g, p13[:, 0], gg); np.add.at(g, p13[:, 1], -gg)
            d = Xc[allpairs[:, 0]]-Xc[allpairs[:, 1]]; r = np.linalg.norm(d, axis=1)
            m = r < 2.7; dd = np.where(m, r-2.7, 0.0); E += (dd**2).sum()
            gg = np.where(m, 2*dd/np.maximum(r, 1e-6), 0.0)[:, None]*d
            np.add.at(g, allpairs[:, 0], gg); np.add.at(g, allpairs[:, 1], -gg)
            return E, g.ravel()
        return lbfgs(energy, X0.ravel()).reshape(-1, 3)
    X = embed_from(cent*R)
    D = np.linalg.norm(X[:, None]-X[None], axis=2)
    for i in range(nvert):
        close = {j for j in range(nvert) if j != i and D[i, j] < 1.8}
        assert close == set(nb[i]), (name, i)
    X -= X.mean(axis=0)
    write(name, ["C"]*nvert, X, bonds)

# ---- nanotubes ----
from ase.build import nanotube
def tube(n, m, cells, name, cap=True, defect=False):
    a = nanotube(n, m, length=cells, bond=1.42)
    X = a.get_positions()
    nC = len(X)
    D = np.linalg.norm(X[:, None]-X[None], axis=2)
    keep = list(range(nC))
    if defect:
        mid = np.argmin(np.abs(X[:, 2]-X[:, 2].mean()) + np.linalg.norm(X[:, :2], axis=1)*0)
        keep.remove(mid)
    X = X[keep]; nC = len(X)
    D = np.linalg.norm(X[:, None]-X[None], axis=2)
    bonds = [(i, j) for i in range(nC) for j in range(i+1, nC) if D[i, j] < 1.7]
    nb = [[] for _ in range(nC)]
    for a_, b_ in bonds: nb[a_].append(b_); nb[b_].append(a_)
    symbols = ["C"]*nC; pos = list(X); allb = list(bonds)
    if cap:
        for i in range(nC):
            if len(nb[i]) == 2:
                d = X[i]-X[nb[i]].mean(axis=0); d /= np.linalg.norm(d)
                pos.append(X[i]+1.09*d); symbols.append("H"); allb.append((i, len(pos)-1))
            else:
                assert len(nb[i]) == 3 or not cap, (name, i, len(nb[i]))
    write(name, symbols, np.array(pos), allb)
tube(5, 5, 14, "nanotube-5-5-hcap")
tube(10, 0, 8, "nanotube-10-0-hcap")
tube(10, 0, 22, "nanotube-10-0-long-hcap")
tube(10, 10, 20, "nanotube-10-10-hcap")
tube(5, 5, 14, "nanotube-5-5-open", cap=False)
tube(10, 0, 8, "nanotube-10-0-open", cap=False)
tube(5, 5, 14, "nanotube-5-5-defect-hcap", defect=True)
