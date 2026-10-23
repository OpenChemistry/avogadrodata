import os, random, sys, urllib.request
from concurrent.futures import ThreadPoolExecutor
from rdkit import Chem
SEED = 20261022
out = sys.argv[1]; os.makedirs(out, exist_ok=True)
curated = """ATP ADP AMP ANP GTP GDP GMP CTP CDP CMP UTP UDP UMP TTP DTP DGT DCT DAT
A G C U T DA DG DC DT ADE GUN CYT URA THM 5MC PSU I
NAD NAP NDP NAI FAD FMN SAM SAH PLP TPP COA ACO HEM CLA BCL HEC HEA PHO
STI IRE ERL LAP 1N1 SUN BAX GIN
CAF NCT FOL MTX BTN RET B12 PQQ F43 SRM
ALA ARG ASP GLU HIS LYS TYR TRP PHE PRO CYS MET
GLC NAG MAN GAL FUC BGC
FES SF4 CU1 HEA CHL BCB 
BNZ NAP IMD TRZ PYR PYZ IND BEN INS TLA
LOV SIM ATV RSV PRA
""".split()
curated = list(dict.fromkeys(curated))

# seeded random sample of the rest, weighted towards (hetero)aromatic ring systems
rows = []
for line in open("ccd/smi.smi"):
    p = line.rstrip("\n").split("\t")
    if len(p) < 2: continue
    rows.append((p[1], p[0]))
rng = random.Random(SEED)
rng.shuffle(rows)
hetero, other = [], []
for cid, smi in rows:
    if cid in curated: continue
    m = Chem.MolFromSmiles(smi)
    if m is None: continue
    n = m.GetNumHeavyAtoms()
    if n < 6 or n > 60: continue
    ri = m.GetRingInfo()
    if ri.NumRings() == 0:
        other.append(cid); continue
    arom_het = any(a.GetIsAromatic() and a.GetSymbol() != "C" for a in m.GetAtoms())
    fused = any(ri.NumAtomRings(a.GetIdx()) > 1 for a in m.GetAtoms())
    (hetero if (arom_het or fused) else other).append(cid)
    if len(hetero) >= 130 and len(other) >= 60: break
picked = hetero[:110] + other[:40]
ids = curated + picked

def fetch(cid):
    try:
        d = urllib.request.urlopen(f"https://files.rcsb.org/ligands/download/{cid}_ideal.sdf", timeout=30).read().decode()
        return cid, d
    except Exception as e:
        return cid, None
kept = []; dropped = []
with ThreadPoolExecutor(8) as ex:
    for cid, d in ex.map(fetch, ids):
        if d is None: dropped.append((cid, "no download")); continue
        lines = d.splitlines()
        try:
            nat = int(lines[3][0:3]); nb = int(lines[3][3:6])
        except Exception:
            dropped.append((cid, "bad header")); continue
        if "V2000" not in lines[3]: dropped.append((cid, "not V2000")); continue
        atoms = lines[4:4+nat]
        heavy = sum(1 for a in atoms if a[31:34].strip() != "H")
        zs = [float(a[20:30]) for a in atoms]
        if heavy > 200 or heavy < 2: dropped.append((cid, f"heavy atoms {heavy}")); continue
        if all(abs(z) < 1e-6 for z in zs): dropped.append((cid, "no z")); continue
        bonds = [int(b[6:9]) for b in lines[4+nat:4+nat+nb]]
        if any(o not in (1, 2, 3) for o in bonds): dropped.append((cid, "aromatic/odd bond type")); continue
        open(os.path.join(out, cid + ".sdf"), "w").write("\n".join(lines[:lines.index("$$$$")] if "$$$$" in lines else lines) + "\n")
        kept.append(cid)
open(os.path.join(out, "ids.txt"), "w").write("\n".join(kept) + "\n")
print(len(ids), "requested,", len(kept), "kept; dropped:", dropped)
print("random picks:", len(picked))
