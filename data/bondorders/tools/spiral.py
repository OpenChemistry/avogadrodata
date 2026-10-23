"""Fullerene graphs from face spirals, and a spring embedding."""
import itertools, sys, random
import numpy as np

def lbfgs(fun, x0, maxiter=5000, m=12, gtol=1e-5):
    x = x0.copy(); f, g = fun(x)
    S, Y = [], []
    for it in range(maxiter):
        if np.abs(g).max() < gtol: break
        q = g.copy(); alphas = []
        for s_, y_ in reversed(list(zip(S, Y))):
            rho = 1.0/np.dot(y_, s_); a = rho*np.dot(s_, q); alphas.append((a, rho, s_, y_)); q -= a*y_
        if S:
            gamma = np.dot(S[-1], Y[-1])/np.dot(Y[-1], Y[-1]); q *= gamma
        else:
            q *= 0.01
        for a, rho, s_, y_ in reversed(alphas):
            b = rho*np.dot(y_, q); q += s_*(a-b)
        d = -q
        if np.dot(d, g) >= 0: d = -g; S, Y = [], []
        t = 1.0; fx = f
        while True:
            xn = x + t*d; fn, gn = fun(xn)
            if fn <= f + 1e-4*t*np.dot(g, d) or t < 1e-12: break
            t *= 0.5
        s_ = xn - x; y_ = gn - g
        if np.dot(s_, y_) > 1e-12:
            S.append(s_); Y.append(y_)
            if len(S) > m: S.pop(0); Y.pop(0)
        x, f, g = xn, fn, gn
    return x

def spiral_dual(degs):
    """Build the dual triangulation from a face-degree spiral.
    Returns (adjacency sets, faces) or None if the spiral is invalid."""
    F = len(degs)
    adj = [set() for _ in range(F)]
    faces = []
    count = [0]*F
    def link(a, b):
        if b in adj[a]: return False
        adj[a].add(b); adj[b].add(a); count[a] += 1; count[b] += 1
        return True
    link(0, 1)
    ring = [0, 1]          # cyclic boundary, oldest-open first
    for v in range(2, F):
        # s = first boundary vertex with open valence; e = newest = last of ring
        s_pos = None
        for i, u in enumerate(ring):
            if degs[u] - count[u] > 0 and u != ring[-1] or (len(ring) == 2 and i == 0):
                s_pos = i; break
        if s_pos is None: return None
        # rotate so that s is first
        ring = ring[s_pos:] + ring[:s_pos]
        s = ring[0]; e = ring[-1]
        if s == e: return None
        if not link(v, e): return None
        if not link(v, s): return None
        faces.append((v, e, s))
        newring = [v] + ring[:]       # v sits between e (last) and s (first)
        # cyclic: ..., e, v, s, ...  -> rotate so list is [s, ..., e, v]
        ring = ring + [v]
        # close triangles at the s side
        while len(ring) > 3 and degs[ring[0]] - count[ring[0]] == 0:
            s0 = ring[0]; s1 = ring[1]
            if s1 == v: break
            if not link(v, s1): return None
            faces.append((v, s0, s1))
            ring = ring[1:]
        for u in range(F):
            if count[u] > degs[u]: return None
    if any(count[u] != degs[u] for u in range(F)):
        return None
    return adj, faces

def fullerene_graph(dual):
    adj, faces = dual
    F = len(adj)
    tri_index = {}
    for i, f in enumerate(faces):
        tri_index[frozenset(f)] = i
    edges = {}
    for i, f in enumerate(faces):
        for a, b in itertools.combinations(f, 2):
            edges.setdefault(frozenset((a, b)), []).append(i)
    bonds = set()
    for k, ts in edges.items():
        if len(ts) == 2:
            bonds.add(tuple(sorted(ts)))
    return len(faces), sorted(bonds), faces

def search_ipr(n, seed=1, limit=None):
    F = n//2 + 2
    # DFS over pentagon positions with spiral validity checked at the end of each prefix
    best = None
    def rec(pent, nextmin):
        nonlocal best
        if best is not None: return
        if len(pent) == 12:
            degs = [5 if (i+1) in pent else 6 for i in range(F)]
            d = spiral_dual(degs)
            if d is None: return
            nv, bonds, faces = fullerene_graph(d)
            if nv != n or len(bonds) != 3*n//2: return
            deg = [0]*nv
            for a, b in bonds: deg[a]+=1; deg[b]+=1
            if set(deg) != {3}: return
            # IPR: pentagon faces (dual deg 5) must not be adjacent
            adj = d[0]
            for p in pent:
                if any((q-0) in pent for q in [x+1 for x in adj[p-1]]): return
            best = (d, sorted(pent))
            return
        for i in range(nextmin, F+1):
            if (12 - len(pent)) > (F - i + 1): break
            rec(pent + [i], i+1)
    rec([], 1)
    return best

def embed(nv, bonds, target, seed=0, tries=40):
    """Spring embedding on a sphere. target: dict bond->length."""
    nb = [[] for _ in range(nv)]
    for a, b in bonds:
        nb[a].append(b); nb[b].append(a)
    pairs13 = set()
    for v in range(nv):
        for a, b in itertools.combinations(nb[v], 2):
            pairs13.add((min(a, b), max(a, b)))
    bset = set(bonds)
    b_idx = np.array(bonds)
    bl = np.array([target[b] for b in bonds])
    p13 = np.array(sorted(pairs13))
    allpairs = np.array([(i, j) for i in range(nv) for j in range(i+1, nv)
                         if (i, j) not in bset and (i, j) not in pairs13])
    rng = np.random.default_rng(seed)
    R = 1.42*np.sqrt(nv/ (4*np.pi) * 2.6)  # rough sphere radius
    def energy(x):
        X = x.reshape(-1, 3)
        g = np.zeros_like(X)
        E = 0.0
        d = X[b_idx[:, 0]] - X[b_idx[:, 1]]
        r = np.linalg.norm(d, axis=1)
        dd = (r - bl)
        E += (dd**2).sum() * 10
        gg = (20*dd/r)[:, None]*d
        np.add.at(g, b_idx[:, 0], gg); np.add.at(g, b_idx[:, 1], -gg)
        d = X[p13[:, 0]] - X[p13[:, 1]]
        r = np.linalg.norm(d, axis=1)
        dd = r - 2.40
        E += (dd**2).sum() * 2
        gg = (4*dd/r)[:, None]*d
        np.add.at(g, p13[:, 0], gg); np.add.at(g, p13[:, 1], -gg)
        if len(allpairs):
            d = X[allpairs[:, 0]] - X[allpairs[:, 1]]
            r = np.linalg.norm(d, axis=1)
            m = r < 2.7
            dd = np.where(m, r - 2.7, 0.0)
            E += (dd**2).sum()
            gg = np.where(m, 2*dd/np.maximum(r, 1e-6), 0.0)[:, None]*d
            np.add.at(g, allpairs[:, 0], gg); np.add.at(g, allpairs[:, 1], -gg)
        return E, g.ravel()
    for t in range(tries):
        X0 = rng.normal(size=(nv, 3)); X0 /= np.linalg.norm(X0, axis=1)[:, None]; X0 *= R
        X = lbfgs(energy, X0.ravel()).reshape(-1, 3)
        # validity: same graph from distances
        D = np.linalg.norm(X[:, None]-X[None], axis=2)
        ok = True
        for i in range(nv):
            close = {j for j in range(nv) if j != i and D[i, j] < 1.8}
            if close != set(nb[i]): ok = False; break
        if ok:
            return X - X.mean(axis=0)
    return None

class Spiral:
    """Incremental spiral builder; degrees are supplied one face at a time."""
    def __init__(self, F):
        self.F = F; self.degs = []; self.adj = []; self.count = []; self.faces = []
        self.ring = []
    def copy(self):
        o = Spiral(self.F); o.degs = self.degs[:]; o.adj = [set(a) for a in self.adj]
        o.count = self.count[:]; o.faces = self.faces[:]; o.ring = self.ring[:]
        return o
    def link(self, a, b):
        if b in self.adj[a]: return False
        self.adj[a].add(b); self.adj[b].add(a); self.count[a] += 1; self.count[b] += 1
        return self.count[a] <= self.degs[a] and self.count[b] <= self.degs[b]
    def add(self, d):
        v = len(self.degs); self.degs.append(d); self.adj.append(set()); self.count.append(0)
        if v == 0: return True
        if v == 1:
            self.ring = [0, 1]; return self.link(0, 1)
        ring = self.ring
        # first boundary vertex with open valence (excluding the newest)
        s_pos = None
        for i, u in enumerate(ring[:-1]):
            if self.degs[u] - self.count[u] > 0:
                s_pos = i; break
        if s_pos is None: return False
        ring = ring[s_pos:] + ring[:s_pos]
        s = ring[0]; e = ring[-1]
        if not (self.link(v, e) and self.link(v, s)): return False
        self.faces.append((v, e, s))
        ring = ring + [v]
        while len(ring) > 3 and self.degs[ring[0]] - self.count[ring[0]] == 0:
            s0, s1 = ring[0], ring[1]
            if s1 == v: break
            if not self.link(v, s1): return False
            self.faces.append((v, s0, s1))
            ring = ring[1:]
        self.ring = ring
        return True
    def done(self):
        return all(self.count[u] == self.degs[u] for u in range(self.F))

def ipr_ok(sp):
    for u in range(len(sp.degs)):
        if sp.degs[u] == 5:
            for w in sp.adj[u]:
                if sp.degs[w] == 5: return False
    return True

def search_ipr2(n, ipr=True):
    F = n//2 + 2
    result = []
    def rec(sp, pents):
        if result: return
        k = len(sp.degs)
        if k == F:
            if pents == 12 and sp.done():
                result.append(sp)
            return
        for d in (5, 6):
            p = pents + (d == 5)
            if p > 12 or (12 - p) > (F - k - 1): continue
            nsp = sp.copy()
            if not nsp.add(d): continue
            if ipr and not ipr_ok(nsp): continue
            rec(nsp, p)
    rec(Spiral(F), 0)
    if not result: return None
    sp = result[0]
    faces = sp.faces + ([tuple(sp.ring)] if len(sp.ring) == 3 else [])
    return (sp.adj, faces), [i+1 for i, d in enumerate(sp.degs) if d == 5]
