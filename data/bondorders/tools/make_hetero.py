import os, sys
from rdkit import Chem
from rdkit.Chem import AllChem
SEED = 0xF00D
out = sys.argv[1]; os.makedirs(out, exist_ok=True)
mols = {
 "pyrrole": "c1cc[nH]c1", "imidazole": "c1c[nH]cn1", "pyrazole": "c1cn[nH]c1",
 "123-triazole-1H": "c1c[nH]nn1", "124-triazole-1H": "[nH]1ncnc1",
 "tetrazole-1H": "[nH]1cnnn1", "tetrazole-2H": "c1nn[nH]n1",
 "indole": "c1ccc2[nH]ccc2c1", "indazole-1H": "c1ccc2[nH]ncc2c1",
 "benzimidazole": "c1ccc2[nH]cnc2c1", "carbazole": "c1ccc2c(c1)[nH]c1ccccc12",
 "purine-7H": "c1ncc2[nH]cnc2n1", "purine-9H": "c1ncc2nc[nH]c2n1",
 "guanine": "Nc1nc2[nH]cnc2c(=O)[nH]1", "cytosine": "Nc1cc[nH]c(=O)n1",
 "uracil": "O=c1cc[nH]c(=O)[nH]1", "thymine": "Cc1c[nH]c(=O)[nH]c1=O",
 "pyrimidin-4-one": "O=c1ccnc[nH]1", "pteridine": "c1cnc2ncncc2n1",
 "phenazine": "c1ccc2nc3ccccc3nc2c1", "acridine": "c1ccc2nc3ccccc3cc2c1",
 "beta-carboline": "c1ccc2c(c1)[nH]c1cnccc12",
 "22-bipyridine": "c1ccc(nc1)-c1ccccn1",
 "110-phenanthroline": "c1cnc2c(c1)ccc1cccnc12",
 "lumiflavin": "Cc1cc2c(cc1C)N(C)C1=NC(=O)NC(=O)C1=N2",
 "adenine-9H": "Nc1ncnc2[nH]cnc12",
}
heavy = ["imidazole", "purine-9H", "adenine-9H", "indole"]
bad = []
for name, smi in mols.items():
    m = Chem.MolFromSmiles(smi)
    if m is None: bad.append((name, "smiles")); continue
    m = Chem.AddHs(m)
    ps = AllChem.ETKDGv3(); ps.randomSeed = SEED
    if AllChem.EmbedMolecule(m, ps) != 0: bad.append((name, "embed")); continue
    try: AllChem.MMFFOptimizeMolecule(m, maxIters=2000)
    except Exception: pass
    Chem.Kekulize(m, clearAromaticFlags=True)
    m.SetProp("_Name", name)
    Chem.MolToMolFile(m, os.path.join(out, name + ".mol"), kekulize=True)
    if name in heavy:
        h = Chem.RemoveHs(m, sanitize=False)
        h.SetProp("_Name", name + "-heavy")
        Chem.MolToMolFile(h, os.path.join(out, name + "-heavy.mol"), kekulize=True)
print(len(os.listdir(out)), "files; skipped:", bad)
