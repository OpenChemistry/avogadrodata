#!/usr/bin/env python3
"""Regenerate the reference-points.json files from critic2 output.

Usage: make_reference_points.py INPUTS_DIR [SF4_DIR]

INPUTS_DIR is the folder holding the case folders h2o/ and n2/, each with
points.xyz and one sub-folder per calculation, each with
critic2/<format>/summary.json (schema avogadro-qtaim-summary/1; lengths in
Bohr, everything else in atomic units). A reference-points.json is written
into each of the data folders listed in CASES, which sit next to this
script's folder. Values are copied exactly as critic2 printed them; nothing is
rounded, converted or recomputed, except that the point positions in Angstrom
are read verbatim from points.xyz.

The one exception is xtb-gfn2: its full QTAIM run is excluded (the density is
valence only), so its values are the point values critic2 printed in
critic2-points/output.cro (the "%% point" blocks, "Field value (f)" and so
on). That output has no K (critic2's virial-based kinetic energy density) and
no per-point summary, so those entries carry rho, gradient, hessian,
laplacian and G only.

With SF4_DIR (holding sf4-points.xyz, sf4.cri and sf4.cro from critic2 1.4.306
run as "load sf4.molden xtb" on data/molden/SF4.molden), also writes
data/molden/SF4-reference-points.json the same way (Angstrom positions only).
"""
import json
import os
import re
import sys

# Each entry: (folder, wavefunction format, file name as stored in the data
# directory, file critic2 actually read, critic2 dialect keyword on the
# "load" line or None, where the values come from).
# The dialect is checked against the load line of the recorded input.cri.
H2O_631GD = [
    ("gaussian-hf-6-31gd-6d", "fchk", "h2o.fchk.gz", "h2o.fchk", None),
    ("gaussian-hf-6-31gd-6d", "wfn", "h2o.wfn", "h2o.wfn", None),
    ("gaussian-hf-6-31gd-6d", "wfx", "h2o.wfx", "h2o.wfx", None),
    ("gaussian-hf-6-31gd-5d", "fchk", "h2o.fchk.gz", "h2o.fchk", None),
    ("gaussian-hf-6-31gd-5d", "wfn", "h2o.wfn", "h2o.wfn", None),
    ("gaussian-hf-6-31gd-5d", "wfx", "h2o.wfx", "h2o.wfx", None),
    ("orca-hf-6-31gd-5d", "molden", "h2o.molden", "h2o.molden", "orca"),
    ("orca-hf-6-31gd-5d", "wfn", "h2o.wfn", "h2o.wfn", None),
    ("orca-hf-6-31gd-5d", "wfx", "h2o.wfx", "h2o.wfx", None),
    # critic2's "pyscf" dialect reads Psi4 1.11 Molden files correctly; its
    # default and "psi4" dialects are 1.8e-2 off for Psi4 1.11 Cartesian d.
    ("psi4-hf-6-31gd-6d", "fchk", "h2o.fchk.gz", "h2o.fchk", None),
    ("psi4-hf-6-31gd-6d", "molden", "h2o.molden", "h2o.molden", "pyscf"),
    ("psi4-hf-6-31gd-5d", "fchk", "h2o.fchk.gz", "h2o.fchk", None),
    ("psi4-hf-6-31gd-5d", "molden", "h2o.molden", "h2o.molden", "pyscf"),
    ("pyscf-hf-6-31gd-6d", "molden", "h2o.molden", "h2o.molden", "pyscf"),
    ("pyscf-hf-6-31gd-5d", "molden", "h2o.molden", "h2o.molden", "pyscf"),
    # no usable QTAIM summary: values from critic2-points/output.cro
    ("xtb-gfn2", "points", "h2o.molden", "molden.input", "xtb"),
]
H2O_CCPVTZ = [
    ("orca-hf-cc-pvtz", "molden", "h2o.molden", "h2o.molden", "orca"),
]
N2_CCPVQZ = [
    ("orca-b3lyp-cc-pvqz", "molden", "n2.molden", "n2.molden", "orca"),
]

CASES = [
    # (data folder, inputs case folder, label, point categories, entries)
    ("h2o-hf-631gd", "h2o", "H2O",
     ["bond-midpoint"] * 2 + ["valence"] * 14 + ["tail"] * 4, H2O_631GD),
    ("h2o-hf-ccpvtz", "h2o", "H2O",
     ["bond-midpoint"] * 2 + ["valence"] * 14 + ["tail"] * 4, H2O_CCPVTZ),
    # n2/points.xyz: "1 midpoint, 15 valence, 4 tail"
    ("n2-b3lyp-ccpvqz", "n2", "N2",
     ["bond-midpoint"] * 1 + ["valence"] * 15 + ["tail"] * 4, N2_CCPVQZ),
]

NUMBER = r"[-+]?\d+\.\d+E[-+]\d+"


def read_points_xyz(path):
    with open(path) as handle:
        lines = handle.read().splitlines()
    count = int(lines[0])
    points = []
    for line in lines[2:2 + count]:
        fields = line.split()
        points.append([float(v) for v in fields[1:4]])
    return points


def load_line(path):
    """Return (file, dialect or None) of the first critic2 "load" command."""
    with open(path) as handle:
        for line in handle:
            fields = line.split()
            if fields and fields[0] == "load":
                return fields[1], (fields[2] if len(fields) > 2 else None)
    raise ValueError("no load command in " + path)


def numbers(text):
    return [float(v) for v in re.findall(NUMBER, text)]


def read_cro_points(path, count):
    """Point values from the "%% point" blocks of a critic2 output file."""
    with open(path) as handle:
        text = handle.read()
    blocks = re.split(r"^%% point.*$", text, flags=re.MULTILINE)[1:]
    assert len(blocks) == count, (path, len(blocks))
    points = []
    for block in blocks:
        # critic2 echoes its input lines ("critic2:15> ...") to the same
        # stream, and one sometimes lands in the middle of a point block,
        # ahead of what was being printed; drop the echo prefix.
        block = re.sub(r"^critic2:\d+>\s?", "", block, flags=re.MULTILINE)

        def field(label):
            match = re.search(r"^\s*" + re.escape(label) + r":(.*)$", block,
                              flags=re.MULTILINE)
            return numbers(match.group(1))

        hessian = re.search(r"^\s*Hessian:\n((?:.*\n){3})", block,
                            flags=re.MULTILINE)
        rows = [numbers(row) for row in hessian.group(1).splitlines()]
        points.append({
            "rho": field("Field value (f)")[0],
            "gradient": field("Gradient (grad f)"),
            "hessian": rows,
            "laplacian": field("Laplacian (del2 f)")[0],
            "G": field("Kinetic energy density (G,tau)")[0],
        })
    return points


def build(source, case):
    data_folder, case_folder, label, categories, entries = case
    case_dir = os.path.join(source, case_folder)
    angstrom = read_points_xyz(os.path.join(case_dir, "points.xyz"))
    assert len(angstrom) == len(categories)

    result = {
        "schema": "avogadro-qtaim-reference-points/1",
        "description": "Electron density and derivatives from critic2 at %d "
                       "fixed points of %s, one entry per wavefunction file."
                       % (len(angstrom), label),
        "units": {
            "length": "Bohr unless the key says angstrom",
            "density": "e/Bohr^3",
            "gradient": "e/Bohr^4",
            "hessian": "e/Bohr^5",
            "laplacian": "e/Bohr^5",
            "G": "Hartree/Bohr^3, as reported by critic2",
            "K": "Hartree/Bohr^3, as reported by critic2",
        },
        "points": [
            {"index": i, "category": categories[i],
             "angstrom": angstrom[i], "bohr": None}
            for i in range(len(angstrom))
        ],
        "provenance": {},
        "wavefunctions": {},
    }

    for folder, fmt, stored, read, dialect in entries:
        if fmt == "points":
            directory = os.path.join(case_dir, folder, "critic2-points")
            cro = os.path.join(directory, "output.cro")
            cri = os.path.join(directory, "input.cri")
            points = read_cro_points(cro, len(angstrom))
            entry = {
                "critic2_input_file": read,
                "critic2_output": "%s/critic2-points/output.cro" % folder,
                "note": "Point values only (no QTAIM run: the xtb density "
                        "is valence only). Parsed from the %% point blocks "
                        "of critic2's output; K is not available.",
                "points": points,
            }
        else:
            directory = os.path.join(case_dir, folder, "critic2", fmt)
            path = os.path.join(directory, "summary.json")
            cri = os.path.join(directory, "input.cri")
            with open(path) as handle:
                summary = json.load(handle)
            meta = summary["metadata"]
            assert summary["schema"] == "avogadro-qtaim-summary/1"
            assert meta["tool"] == "critic2"
            assert not summary.get("errors"), (path, summary.get("errors"))
            assert len(summary["points"]) == len(angstrom)

            prov = result["provenance"]
            for key in ("tool", "tool_version", "bohr_angstrom"):
                assert prov.setdefault(key, meta[key]) == meta[key], (path, key)

            for i, point in enumerate(summary["points"]):
                ref = result["points"][i]
                if ref["bohr"] is None:
                    ref["bohr"] = point["position"]
                assert ref["bohr"] == point["position"], (path, i)

            entry = {
                "critic2_input_file": read,
                "critic2_summary": "%s/critic2/%s/summary.json" % (folder, fmt),
                "electrons": summary["electrons"],
                "points": [
                    {name: point[name]
                     for name in ("rho", "gradient", "hessian", "laplacian",
                                  "G", "K") if name in point}
                    for point in summary["points"]
                ],
            }
        # the dialect the reference was computed with, from its own input
        loaded, loaded_dialect = load_line(cri)
        assert loaded_dialect == dialect, (cri, loaded_dialect, dialect)
        if dialect is not None:
            entry["critic2_dialect"] = dialect
        result["wavefunctions"]["%s/%s" % (folder, stored)] = entry

    # the xtb output only has positions in Angstrom; the Bohr positions come
    # from the summaries of the other entries of the same case
    for point in result["points"]:
        assert point["bohr"] is not None
    return result


def build_sf4(directory):
    angstrom = read_points_xyz(os.path.join(directory, "sf4-points.xyz"))
    cri = os.path.join(directory, "sf4.cri")
    loaded, dialect = load_line(cri)
    assert (loaded, dialect) == ("sf4.molden", "xtb"), (loaded, dialect)
    with open(os.path.join(directory, "sf4.cro")) as handle:
        version = re.search(r"critic2 \(.*\), version (\S+)",
                            handle.read()).group(1)
    return {
        "schema": "avogadro-qtaim-reference-points/1",
        "description": "Electron density and derivatives from critic2 at %d "
                       "fixed points around S of SF4, for data/molden/"
                       "SF4.molden (xtb, Cartesian d)." % len(angstrom),
        "units": {
            "length": "Angstrom",
            "density": "e/Bohr^3",
            "gradient": "e/Bohr^4",
            "hessian": "e/Bohr^5",
            "laplacian": "e/Bohr^5",
            "G": "Hartree/Bohr^3, as reported by critic2",
        },
        "points": [{"index": i, "angstrom": p}
                   for i, p in enumerate(angstrom)],
        "provenance": {"tool": "critic2", "tool_version": version},
        "wavefunctions": {
            "SF4.molden": {
                "critic2_input_file": "sf4.molden",
                "critic2_dialect": "xtb",
                "critic2_output": "point blocks of critic2 output "
                                  "(sf4.cro), not stored",
                "points": read_cro_points(
                    os.path.join(directory, "sf4.cro"), len(angstrom)),
            }
        },
    }


def main():
    source = sys.argv[1]
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here)
    for case in CASES:
        result = build(source, case)
        output = os.path.join(root, case[0], "reference-points.json")
        with open(output, "w") as handle:
            json.dump(result, handle, indent=1)
            handle.write("\n")
    if len(sys.argv) > 2:
        output = os.path.join(os.path.dirname(root), "molden",
                              "SF4-reference-points.json")
        with open(output, "w") as handle:
            json.dump(build_sf4(sys.argv[2]), handle, indent=1)
            handle.write("\n")


if __name__ == "__main__":
    main()
