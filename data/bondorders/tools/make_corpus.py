import os, sys
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

SEED = 0xF00D
out = sys.argv[1]
os.makedirs(out, exist_ok=True)
mols = {
 "benzene": "c1ccccc1", "toluene": "Cc1ccccc1", "phenol": "Oc1ccccc1",
 "aniline": "Nc1ccccc1", "pyridine": "c1ccncc1", "pyrrole": "c1cc[nH]c1",
 "furan": "c1ccoc1", "thiophene": "c1ccsc1", "imidazole": "c1c[nH]cn1",
 "pyrimidine": "c1cncnc1", "naphthalene": "c1ccc2ccccc2c1",
 "anthracene": "c1ccc2cc3ccccc3cc2c1", "phenanthrene": "c1ccc2c(c1)ccc1ccccc12",
 "pyrene": "c1cc2ccc3cccc4ccc(c1)c2c34",
 "tetracene": "c1ccc2cc3cc4ccccc4cc3cc2c1",
 "coronene": "c1cc2ccc3ccc4ccc5ccc6ccc1c1c2c3c4c5c61",
 "porphine": "C1=CC2=NC1=CC1=CC=C(N1)C=C1C=CC(=N1)C=C1C=CC(=C2)N1",
 "indole": "c1ccc2[nH]ccc2c1", "quinoline": "c1ccc2ncccc2c1",
 "adenine": "Nc1ncnc2[nH]cnc12",
 "caffeine": "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
 "2-pyridone": "O=c1cccc[nH]1", "4-pyridone": "O=c1cc[nH]cc1",
 "azobenzene": "c1ccc(cc1)N=Nc1ccccc1",
 "stilbene": "C(=Cc1ccccc1)c1ccccc1", "biphenyl": "c1ccc(cc1)-c1ccccc1",
 "acetone": "CC(C)=O", "formaldehyde": "C=O", "acetic-acid": "CC(O)=O",
 "acetamide": "CC(N)=O", "benzamide": "NC(=O)c1ccccc1",
 "acetonitrile": "CC#N", "benzonitrile": "N#Cc1ccccc1",
 "acetylene": "C#C", "ethylene": "C=C", "butadiene": "C=CC=C",
 "allene": "C=C=C", "acrolein": "C=CC=O", "ethyl-acetate": "CCOC(C)=O",
 "urea": "NC(N)=O", "cyclohexane": "C1CCCCC1", "adamantane": "C1C2CC3CC1CC(C2)C3",
 "ethanol": "CCO", "cyclohexene": "C1=CCCCC1", "norbornane": "C1CC2CCC1C2",
 "aspirin": "CC(=O)Oc1ccccc1C(O)=O", "paracetamol": "CC(=O)Nc1ccc(O)cc1",
 "glycine": "NCC(O)=O", "isoprene": "C=CC(C)=C", "dmf": "CN(C)C=O",
 "azulene": "c1ccc2cccc2cc1", "chlorobenzene": "Clc1ccccc1",
 "trifluorotoluene": "FC(F)(F)c1ccccc1",
 "nicotine": "CN1CCCC1c1cccnc1",
 "ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(O)=O",
 "hexatriene": "C=CC=CC=C", "cyclopentadiene": "C1=CC=CC1",
 "benzoquinone": "O=C1C=CC(=O)C=C1", "pyrazine": "c1cnccn1",
 "oxazole": "c1cocn1", "thiazole": "c1cscn1",
 "acetophenone": "CC(=O)c1ccccc1", "benzaldehyde": "O=Cc1ccccc1",
 "diphenylacetylene": "c1ccc(cc1)C#Cc1ccccc1",
 "carbon-disulfide": "S=C=S", "hcn": "C#N", "ketene": "C=C=O",
}
bad = []
for name, smi in mols.items():
    m = Chem.MolFromSmiles(smi)
    if m is None:
        bad.append((name, "smiles")); continue
    m = Chem.AddHs(m)
    ps = AllChem.ETKDGv3(); ps.randomSeed = SEED
    if AllChem.EmbedMolecule(m, ps) != 0:
        bad.append((name, "embed")); continue
    try:
        AllChem.MMFFOptimizeMolecule(m, maxIters=2000)
    except Exception:
        pass
    Chem.Kekulize(m, clearAromaticFlags=True)
    m.SetProp("_Name", name)
    Chem.MolToMolFile(m, os.path.join(out, name + ".mol"), kekulize=True)
print(len(mols) - len(bad), "written; skipped:", bad)
print("RDKit", rdBase.rdkitVersion)
