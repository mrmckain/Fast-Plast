#!/usr/bin/env python3
"""
repair_taxonomy.py

Re-fetch taxonomy for plastome accessions whose ranks are missing or garbled in
the Fast-Plast metadata table, using NCBI Entrez (accession -> TaxId -> lineage),
and write a corrected metadata TSV.

This is the companion to normalize_plastome_db.pl. When the original taxonomy
pull failed for some records, their rows land with NA/blank ranks; this repairs
exactly those rows (or an explicit accession list) without touching the rest.

Fills the ranks Fast-Plast uses: genus, species, tribe, subfamily, family, order.
Ranks that NCBI genuinely does not assign for a lineage are left as 'NA' -- this
does not invent taxonomy, it only recovers what NCBI actually has.

Requires Biopython (pip install biopython) and network access to NCBI E-utilities.
NCBI requires an email; an API key (free) raises the rate limit from 3 to 10/sec.

Usage:
  repair_taxonomy.py --meta GenBank_Plastomes.metadata.tsv \\
                     --out  GenBank_Plastomes.metadata.fixed.tsv \\
                     --email you@example.edu [--api-key KEY] \\
                     [--accessions to_fix.txt]

Default: repairs every row whose 'order' is NA/blank.
With --accessions: repairs exactly the accessions listed (one per line).
"""

import argparse
import sys
import time

COLUMNS = ["accession", "genus", "species", "tribe", "subfamily", "family", "order", "source"]
WANT_RANKS = {"genus", "tribe", "subfamily", "family", "order"}


def ranks_from_lineage(organism_name, lineage_ex):
    """Pure function: build the Fast-Plast rank fields from an NCBI taxonomy record.

    organism_name: e.g. "Arabidopsis thaliana"
    lineage_ex:    list of dicts each with 'ScientificName' and 'Rank'
                   (Biopython's LineageEx structure)
    Returns dict with genus/species/tribe/subfamily/family/order, 'NA' where absent.
    """
    out = {r: "NA" for r in ("genus", "species", "tribe", "subfamily", "family", "order")}

    for node in lineage_ex:
        rank = str(node.get("Rank", "")).lower()
        name = str(node.get("ScientificName", "")).strip()
        if rank in WANT_RANKS and name:
            out[rank] = name

    # species epithet from the binomial (schema stores the epithet, not "Genus species")
    toks = (organism_name or "").split()
    if len(toks) >= 2:
        if out["genus"] == "NA":
            out["genus"] = toks[0]
        out["species"] = toks[1]
    elif len(toks) == 1 and out["genus"] == "NA":
        out["genus"] = toks[0]

    return out


def load_meta(path):
    rows = []  # list of dicts, preserving order
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if not line.strip():
                continue
            vals = line.rstrip("\n").split("\t")
            vals += [""] * (len(header) - len(vals))
            rows.append(dict(zip(header, vals)))
    return header, rows


def needs_repair(row):
    o = row.get("order", "").strip()
    return o == "" or o.upper() == "NA"


def chunked(seq, n):
    for i in range(0, len(seq), n):
        yield seq[i:i + n]


def main():
    ap = argparse.ArgumentParser(description="Repair Fast-Plast taxonomy metadata via NCBI.")
    ap.add_argument("--meta", required=True, help="input metadata TSV")
    ap.add_argument("--out", required=True, help="output (corrected) metadata TSV")
    ap.add_argument("--email", required=True, help="your email (required by NCBI)")
    ap.add_argument("--api-key", default=None, help="NCBI API key (optional, higher rate limit)")
    ap.add_argument("--accessions", default=None,
                    help="file of accessions to repair (default: rows with order=NA)")
    ap.add_argument("--batch", type=int, default=50, help="ids per Entrez request")
    args = ap.parse_args()

    try:
        from Bio import Entrez
    except ImportError:
        sys.exit("ERROR: Biopython not found. Install with: pip install biopython")

    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key
    delay = 0.11 if args.api_key else 0.34  # stay under 10/s or 3/s

    header, rows = load_meta(args.meta)
    by_acc = {r["accession"]: r for r in rows}

    if args.accessions:
        with open(args.accessions) as fh:
            targets = [ln.strip() for ln in fh if ln.strip()]
        targets = [a for a in targets if a in by_acc]
    else:
        targets = [r["accession"] for r in rows if needs_repair(r)]

    if not targets:
        sys.stderr.write("Nothing to repair; writing metadata unchanged.\n")
    else:
        sys.stderr.write("Repairing %d accession(s) via NCBI...\n" % len(targets))

    # accession normalizer: strip version + uppercase, so "NC_000932.1",
    # "NC_000932", and "nc_000932.1" all compare equal. esummary sometimes
    # returns the bare Caption (no version) while our rows carry the version;
    # matching on the normalized form avoids silently missing every row.
    def norm(a):
        return str(a).split(".")[0].strip().upper()

    # 1) accession -> TaxId  (esummary on nuccore); keyed by normalized accession
    acc2taxid = {}
    for batch in chunked(targets, args.batch):
        h = Entrez.esummary(db="nuccore", id=",".join(batch), post=True)
        recs = Entrez.read(h); h.close()
        for rec in recs:
            acc = str(rec.get("AccessionVersion") or rec.get("Caption") or "")
            taxid = str(rec.get("TaxId", "") or "")
            if acc and taxid and taxid != "0":
                acc2taxid[norm(acc)] = taxid
        time.sleep(delay)

    resolved = sum(1 for a in targets if norm(a) in acc2taxid)
    sys.stderr.write("  resolved %d/%d accessions to TaxIds\n" % (resolved, len(targets)))
    if resolved == 0:
        # expose the actual docsum keys so a remaining mismatch is diagnosable
        try:
            h = Entrez.esummary(db="nuccore", id=targets[0]); first = Entrez.read(h)[0]; h.close()
            sys.stderr.write("  DEBUG docsum keys for %s: %s\n"
                             % (targets[0], ", ".join(sorted(map(str, first.keys())))))
        except Exception as e:
            sys.stderr.write("  DEBUG could not introspect docsum: %s\n" % e)

    # 2) TaxId -> lineage  (efetch on taxonomy)
    taxids = sorted(set(acc2taxid.values()))
    taxid2ranks = {}
    for batch in chunked(taxids, args.batch):
        h = Entrez.efetch(db="taxonomy", id=",".join(batch), post=True)
        recs = Entrez.read(h); h.close()
        for rec in recs:
            tid = str(rec.get("TaxId", ""))
            organism = str(rec.get("ScientificName", ""))
            lineage_ex = rec.get("LineageEx", [])
            # include the organism's own node in case its rank (e.g. species) is useful
            self_node = {"ScientificName": organism, "Rank": str(rec.get("Rank", ""))}
            taxid2ranks[tid] = ranks_from_lineage(organism, list(lineage_ex) + [self_node])
        time.sleep(delay)
    sys.stderr.write("  fetched lineage for %d/%d TaxIds\n" % (len(taxid2ranks), len(taxids)))

    # 3) patch rows (look up via normalized accession)
    fixed = 0
    for acc in targets:
        tid = acc2taxid.get(norm(acc))
        if not tid or tid not in taxid2ranks:
            continue
        r = by_acc[acc]
        newvals = taxid2ranks[tid]
        changed = False
        for k in ("genus", "species", "tribe", "subfamily", "family", "order"):
            if newvals.get(k, "NA") != "NA" and r.get(k, "") != newvals[k]:
                r[k] = newvals[k]
                changed = True
        if changed:
            fixed += 1

    with open(args.out, "w") as out:
        out.write("\t".join(header) + "\n")
        for r in rows:
            out.write("\t".join(r.get(c, "") for c in header) + "\n")

    still_na = sum(1 for a in targets
                   if by_acc[a].get("order", "").strip().upper() in ("", "NA"))
    sys.stderr.write("Done. Rows updated: %d. Still lacking order after repair: %d.\n"
                     % (fixed, still_na))
    if still_na:
        sys.stderr.write("  (those TaxIds genuinely have no assigned order in NCBI, "
                         "or the accession was not resolvable.)\n")


if __name__ == "__main__":
    main()
