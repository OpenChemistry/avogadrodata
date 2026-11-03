#!/usr/bin/env python3
"""Regenerate reference-points.json from critic2 summary.json files.

Usage: make_reference_points.py INPUTS_H2O_DIR [OUTPUT]

INPUTS_H2O_DIR is the folder holding points.xyz and one sub-folder per
calculation, each with critic2/<format>/summary.json (schema
avogadro-qtaim-summary/1; lengths in Bohr, everything else in atomic units).
Values are copied exactly as critic2 printed them; nothing is rounded,
converted or recomputed, except that the point positions in Angstrom are
read verbatim from points.xyz.
"""
import json
import os
import sys

# (folder, wavefunction format, file name as stored in this directory,
#  file critic2 actually read)
ENTRIES = [
    ("gaussian-hf-6-31gd-6d", "fchk", "h2o.fchk.gz", "h2o.fchk"),
    ("gaussian-hf-6-31gd-6d", "wfn", "h2o.wfn", "h2o.wfn"),
    ("gaussian-hf-6-31gd-6d", "wfx", "h2o.wfx", "h2o.wfx"),
    ("gaussian-hf-6-31gd-5d", "fchk", "h2o.fchk.gz", "h2o.fchk"),
    ("gaussian-hf-6-31gd-5d", "wfn", "h2o.wfn", "h2o.wfn"),
    ("gaussian-hf-6-31gd-5d", "wfx", "h2o.wfx", "h2o.wfx"),
    ("orca-hf-6-31gd-5d", "molden", "h2o.molden", "h2o.molden"),
    ("orca-hf-6-31gd-5d", "wfn", "h2o.wfn", "h2o.wfn"),
    ("orca-hf-6-31gd-5d", "wfx", "h2o.wfx", "h2o.wfx"),
]
# points.xyz comment line: "2 midpoint, 14 valence, 4 tail"
CATEGORIES = ["bond-midpoint"] * 2 + ["valence"] * 14 + ["tail"] * 4


def read_points_xyz(path):
    with open(path) as handle:
        lines = handle.read().splitlines()
    count = int(lines[0])
    points = []
    for line in lines[2:2 + count]:
        fields = line.split()
        points.append([float(v) for v in fields[1:4]])
    return points


def main():
    source = sys.argv[1]
    output = sys.argv[2] if len(sys.argv) > 2 else "reference-points.json"
    angstrom = read_points_xyz(os.path.join(source, "points.xyz"))
    assert len(angstrom) == len(CATEGORIES)

    result = {
        "schema": "avogadro-qtaim-reference-points/1",
        "description": "Electron density and derivatives from critic2 at 20 "
                       "fixed points of H2O, one entry per wavefunction file.",
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
            {"index": i, "category": CATEGORIES[i],
             "angstrom": angstrom[i], "bohr": None}
            for i in range(len(angstrom))
        ],
        "provenance": {},
        "wavefunctions": {},
    }

    for folder, fmt, stored, read in ENTRIES:
        path = os.path.join(source, folder, "critic2", fmt, "summary.json")
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

        key = "%s/%s" % (folder, stored)
        result["wavefunctions"][key] = {
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

    with open(output, "w") as handle:
        json.dump(result, handle, indent=1)
        handle.write("\n")


if __name__ == "__main__":
    main()
