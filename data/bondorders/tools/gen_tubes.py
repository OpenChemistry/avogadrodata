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

# ---- nanotubes ----
from math import gcd
def roll(n, m, cells, bond=1.42):
    a = bond*np.sqrt(3)
    a1 = np.array([a, 0.0]); a2 = np.array([a/2, a*np.sqrt(3)/2])
    C = n*a1 + m*a2
    dR = gcd(2*m+n, 2*n+m)
    T = ((2*m+n)//dR)*a1 - ((2*n+m)//dR)*a2
    Lc = np.linalg.norm(C); Lt = np.linalg.norm(T)
    ec = C/Lc; et = T/Lt
    N = 2*(abs(n)+abs(m)) + cells*(abs(2*m+n)+abs(2*n+m))//dR + 6
    I, J = np.meshgrid(np.arange(-N, N), np.arange(-N, N), indexing="ij")
    I = I.ravel(); J = J.ravel()
    pts = []
    for off in (np.zeros(2), (a1+a2)/3):
        r = I[:, None]*a1 + J[:, None]*a2 + off
        u = r@ec/Lc; v = r@et
        k = (u >= -1e-6) & (u < 1-1e-6) & (v >= -1e-6) & (v < cells*Lt-1e-6)
        ang = 2*np.pi*u[k]; R = Lc/(2*np.pi)
        pts.extend(zip(R*np.cos(ang), R*np.sin(ang), v[k]))
    return np.array(pts)

def tube(n, m, cells, name, cap=True, defect=False):
    X = roll(n, m, cells)
    nC = len(X)
    keep = list(range(nC))
    if defect:
        mid = np.argmin(np.abs(X[:, 2]-X[:, 2].mean()) + np.linalg.norm(X[:, :2], axis=1)*0)
        keep.remove(mid)
    X = X[keep]; nC = len(X)
    while True:  # drop atoms that would dangle by a single bond
        D = np.linalg.norm(X[:, None]-X[None], axis=2)
        deg = ((D < 1.7).sum(axis=1)) - 1
        if (deg >= 2).all(): break
        X = X[deg >= 2]; nC = len(X)
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
