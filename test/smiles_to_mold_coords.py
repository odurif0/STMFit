#!/usr/bin/env python3

"""Generate STMFit mold coordinate TSVs from mapped SMILES using RDKit.

This is an optional helper. It is intentionally not a Julia dependency: run it
only in an environment where RDKit is installed. Use atom-map numbers in the
SMILES to identify anchor atoms, for example C1/C2/C4 can be encoded as map
numbers 1/2/4 and passed with ``--anchor-map C1:1,C2:2,C4:4``.

Output columns are consumed by ``test/project_mold_atoms.jl``:
  type, atom, element, x_nm, y_nm, z_nm
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


def parse_anchor_map(text: str) -> dict[int, str]:
    out: dict[int, str] = {}
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        if ":" not in item:
            raise ValueError(f"invalid anchor mapping {item!r}; expected NAME:MAPNUM")
        name, num = item.split(":", 1)
        out[int(num)] = name.strip()
    return out


def load_smiles(args: argparse.Namespace) -> list[tuple[int, str]]:
    pairs: list[tuple[int, str]] = []
    if args.smiles_tsv:
        with open(args.smiles_tsv, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                if not row or row.get("type", "").startswith("#"):
                    continue
                pairs.append((int(row["type"]), row["smiles"].strip()))
    if args.smiles0:
        pairs.append((0, args.smiles0))
    if args.smiles1:
        pairs.append((1, args.smiles1))
    if not pairs:
        raise ValueError("provide --smiles-tsv or both --smiles0/--smiles1")
    types = {typ for typ, _ in pairs}
    missing = {0, 1} - types
    if missing:
        raise ValueError(f"missing SMILES for type(s): {sorted(missing)}")
    return sorted(pairs, key=lambda x: x[0])


def embed_molecule(smiles: str, seed: int, max_iters: int):
    try:
        from rdkit import Chem  # type: ignore[import-not-found]
        from rdkit.Chem import AllChem  # type: ignore[import-not-found]
    except ImportError as exc:
        raise SystemExit(
            "RDKit is not installed. Install it in this Python environment, e.g. "
            "`conda install -c conda-forge rdkit`, then rerun this script."
        ) from exc

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    status = AllChem.EmbedMolecule(mol, params)
    if status != 0:
        raise ValueError(f"RDKit embedding failed for SMILES: {smiles}")
    if AllChem.MMFFHasAllMoleculeParams(mol):
        AllChem.MMFFOptimizeMolecule(mol, maxIters=max_iters)
    else:
        AllChem.UFFOptimizeMolecule(mol, maxIters=max_iters)
    return mol


def atom_name(atom, anchor_by_map: dict[int, str]) -> str:
    mapnum = atom.GetAtomMapNum()
    if mapnum in anchor_by_map:
        return anchor_by_map[mapnum]
    return f"{atom.GetSymbol()}{atom.GetIdx()}"


def write_coords(args: argparse.Namespace) -> None:
    anchor_by_map = parse_anchor_map(args.anchor_map)
    pairs = load_smiles(args)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["type", "atom", "element", "x_nm", "y_nm", "z_nm"])
        for typ, smiles in pairs:
            mol = embed_molecule(smiles, args.seed + typ, args.max_iters)
            conf = mol.GetConformer()
            seen_names: set[str] = set()
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() == 1 and not args.include_hydrogens:
                    continue
                name = atom_name(atom, anchor_by_map)
                if name in seen_names:
                    name = f"{name}_{atom.GetIdx()}"
                seen_names.add(name)
                pos = conf.GetAtomPosition(atom.GetIdx())
                writer.writerow([
                    typ,
                    name,
                    atom.GetSymbol(),
                    f"{pos.x / 10.0:.8g}",
                    f"{pos.y / 10.0:.8g}",
                    f"{pos.z / 10.0:.8g}",
                ])


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Generate chitosan_mold_coords.tsv from mapped SMILES using RDKit."
    )
    p.add_argument("--smiles-tsv", help="TSV with columns: type, smiles")
    p.add_argument("--smiles0", help="Mapped SMILES for type 0 (GlcN)")
    p.add_argument("--smiles1", help="Mapped SMILES for type 1 (GlcNAc)")
    p.add_argument("--out", default="templates/chitosan_mold_coords.tsv")
    p.add_argument(
        "--anchor-map",
        default="C1:1,C2:2,C4:4",
        help="Comma-separated NAME:ATOMMAP mapping for frame anchors [C1:1,C2:2,C4:4]",
    )
    p.add_argument("--include-hydrogens", action="store_true")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--max-iters", type=int, default=500)
    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        write_coords(args)
    except Exception as exc:  # keep CLI failures readable for non-Python users
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(f"Wrote: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
