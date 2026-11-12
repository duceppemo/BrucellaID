#!/usr/bin/env python3
"""Rebuild the list of public read sets the at-scale sweep reports as "Brucella sp.".

The sweep output, the submitter labels and the SKESA assemblies are not shipped (see
validation/at-scale-notes.md), so this takes their locations as arguments and defaults to where
they were kept. It writes one row per run with the sweep's call, the submitter's label and the
assembly's size, and marks a run "clean" when its assembly has at most 150 contigs and is within
0.5 Mb of 3.3 Mb -- a criterion that fits the classical species only: a former Ochrobactrum
genome is 4.4-5.9 Mb and fails it by design, not by quality.

    python3 -I curation/unnamed_lineages/make_candidates.py > candidates_raw.tsv
"""
import argparse
import csv
import os
import sys

A = "/media/8tb/analyses"


def table(path):
    with open(path, newline="") as handle:
        rows = [r for r in csv.reader(handle, delimiter="\t") if r and not r[0].startswith("#")]
    return [dict(zip(rows[0], r)) for r in rows[1:]]


def assembly_stats(fasta):
    lengths = []
    with open(fasta) as handle:
        for line in handle:
            if line.startswith(">"):
                lengths.append(0)
            else:
                lengths[-1] += len(line.strip())
    lengths.sort(reverse=True)
    total, running = sum(lengths), 0
    for length in lengths:
        running += length
        if running >= total / 2:
            return len(lengths), total, length
    return 0, 0, 0


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--extended", default=f"{A}/brucellaid_211_resweep/results_extended_211.tsv")
    p.add_argument("--default", default=f"{A}/brucellaid_211_resweep/results_default_211.tsv")
    p.add_argument("--labels", default=f"{A}/brucellaid_large-scale_validation/phase2/runs.tsv")
    # the three-file runs first: their paired-only assembly, where one exists, is the mis-paired one
    p.add_argument("--assemblies", nargs="+", default=[f"{A}/kmer_rederive/assemblies_unpaired",
                                                       f"{A}/kmer_rederive/assemblies"])
    p.add_argument("--max-contigs", type=int, default=150)
    p.add_argument("--size", type=float, default=3.3e6)
    p.add_argument("--size-window", type=float, default=0.5e6)
    args = p.parse_args()

    default = {r["Sample"]: r for r in table(args.default)}
    label = {r["run"]: r["species"] for r in table(args.labels)}
    out = []
    for r in table(args.extended):
        if r["Species"] != "Brucella sp.":
            continue
        run = r["Sample"]
        fasta = next((f"{d}/{run}.fasta" for d in args.assemblies
                      if os.path.exists(f"{d}/{run}.fasta") and os.path.getsize(f"{d}/{run}.fasta")), "")
        contigs, total, n50 = assembly_stats(fasta) if fasta else ("", "", "")
        clean = bool(fasta) and contigs <= args.max_contigs and abs(total - args.size) <= args.size_window
        d = default.get(run, {})
        out.append({
            "run": run,
            "submitted_as": label.get(run, ""),
            "closest_reference": "_".join(r["ClosestReference"].split("_")[:2]),
            "mash_distance": r["MashDistance"],
            "shared_hashes": r["SharedHashes"],
            "sequence_type": r["SequenceType"],
            "allele_profile": r["AlleleProfile"],
            "default_set_call": d.get("Species", ""),
            "default_set_distance": d.get("MashDistance", ""),
            "contigs": contigs,
            "assembly_bp": total,
            "n50": n50,
            "clean_and_right_sized": "yes" if clean else "no",
            "assembly": fasta,
        })
    out.sort(key=lambda x: (x["closest_reference"], float(x["mash_distance"])))
    writer = csv.DictWriter(sys.stdout, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(out)
    print(f"{len(out)} runs called Brucella sp., "
          f"{sum(o['clean_and_right_sized'] == 'yes' for o in out)} clean and right-sized",
          file=sys.stderr)


if __name__ == "__main__":
    main()
