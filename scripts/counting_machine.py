import os, glob, argparse, sys, json, re

#Counts gene/transcript models and creates a summary table
_SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))
if _SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, _SCRIPTS_DIR)
from buscoPlot import run as _busco_plot


def run_busco_plot(glob_pattern, output_path):
    _busco_plot(glob_pattern, output_path)


BUSCO_FLAT_COLS = ["species", "src", "lineage_used", "lineage_completeness", "eukaryote_completeness"]

# prediction files are named <species>_<taxid>_<source>.<ext>; the source (annotation
# origin) is parsed from the name, so categories only need the directory + extension.
CATEGORIES = {
    "regular": {"pred_dir": "pred",    "pred_ext": ".gff3"},
    "merged":  {"pred_dir": "merged",  "pred_ext": ".gff"},
}

# one consolidated row per species (same contract as the annotator family's
# _metrics.tsv; see isoquant_annotator/derived_metrics.md). geneid is an
# ab-initio coding-gene predictor: it builds no transcriptome and runs no ORF
# calling, so coding_transcripts / transcriptome_transcripts / coding_fraction
# are not measured here and stay NA. Density and isoforms/gene follow from the
# model counts + assembly size.
METRIC_HEADER = [
    "species",
    "gene_count",
    "transcript_count",
    "genome_size_bp",
    "coding_transcripts",
    "transcriptome_transcripts",
    "gene_density_per_mb",
    "transcript_density_per_mb",
    "isoforms_per_gene",
    "coding_fraction",
]


def derived_metrics(genes, transcripts, gsize):
    """(gene_density_per_mb, transcript_density_per_mb, isoforms_per_gene) as
    preformatted strings; NA when a denominator is zero/None."""
    gd = f"{genes / (gsize / 1e6):.2f}" if gsize else "NA"
    td = f"{transcripts / (gsize / 1e6):.2f}" if gsize else "NA"
    ipg = f"{transcripts / genes:.3f}" if genes else "NA"
    return gd, td, ipg


def count_models(gff_path):
    """Count gene and transcript (mRNA) models in a GFF/GFF3 file."""
    genes = 0
    transcripts = 0
    with open(gff_path) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 3:
                continue
            feature = cols[2]
            if feature == "gene":
                genes += 1
            elif feature == "mRNA":
                transcripts += 1
    return genes, transcripts


def genome_size_bp(species_root, sp):
    """Total assembly length (exact; sum of contig lengths, incl. N gaps) for sp.

    Reads the CLEAN_*.fna genome under <species_root>/<sp>/. Prefers the .fai
    index if present; otherwise sums sequence line lengths from the FASTA.
    Returns the size as an int, or None if no genome FASTA is found.
    """
    hits = sorted(glob.glob(os.path.join(species_root, sp, "CLEAN_*.fna")))
    if not hits:
        return None
    fasta = hits[0]
    fai = fasta + ".fai"
    if os.path.isfile(fai):
        total = 0
        with open(fai) as fh:
            for line in fh:
                parts = line.split("\t")
                if len(parts) >= 2:
                    total += int(parts[1])
        return total
    total = 0
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            total += len(line.strip())
    return total


def write_counts(results_dir, summary_dir, name, cfg, species_root=None):
    """Write per-species counts + derived metrics for one evaluation flavour into
    results/summary/<name>/: a counts/ folder with one consolidated <sp>_metrics.tsv
    per species, plus a counts_summary.tsv table (same columns) at the category root.

    Prediction files are named <species>_<source>.<ext> (source = annotation origin).
    The per-species outputs keep one row per species (cross-project _metrics.tsv
    contract); when a species has several sources the first (sorted) one is its
    representative. Per-(species, source) counts are also collected for the local
    summary. Genome size is the total assembly length read from
    <species_root>/<sp>/CLEAN_*.fna; species_root defaults to
    <results_dir>/../training_data/species.

    Returns (rows, unit_rows): rows is the per-species list of metric lists;
    unit_rows is a list of dicts {unit, sp, src, metrics}, one per prediction file.
    """
    if species_root is None:
        species_root = os.path.join(os.path.dirname(os.path.abspath(results_dir)),
                                    "training_data", "species")
    cat_out = os.path.join(summary_dir, name)
    counts_out = os.path.join(cat_out, "counts")
    os.makedirs(counts_out, exist_ok=True)

    pred_dir = os.path.join(results_dir, cfg["pred_dir"])
    ext = cfg["pred_ext"]
    rows = []
    unit_rows = []
    seen_sp = set()
    for gff in sorted(glob.glob(os.path.join(pred_dir, "*" + ext))):
        unit = os.path.basename(gff)[: -len(ext)]  # <specie_name>_<taxid>_<source>
        sp, src = unit.rsplit("_", 1)
        genes, transcripts = count_models(gff)
        gsize = genome_size_bp(species_root, sp)
        gd, td, ipg = derived_metrics(genes, transcripts, gsize)
        metrics = [sp, str(genes), str(transcripts),
                   "NA" if gsize is None else str(gsize),
                   "NA", "NA", gd, td, ipg, "NA"]
        unit_rows.append({"unit": unit, "sp": sp, "src": src, "metrics": metrics})
        # per-species outputs keep one row per species (representative = first source)
        if sp not in seen_sp:
            seen_sp.add(sp)
            with open(os.path.join(counts_out, f"{sp}_metrics.tsv"), "w") as fh:
                fh.write("\t".join(METRIC_HEADER) + "\n")
                fh.write("\t".join(metrics) + "\n")
            rows.append(metrics)

    with open(os.path.join(cat_out, "counts_summary.tsv"), "w") as fh:
        fh.write("\t".join(METRIC_HEADER) + "\n")
        for row in rows:
            fh.write("\t".join(row) + "\n")
    print(f"[{name}] metrics: {len(rows)} species ({len(unit_rows)} units) -> counts/ (+ counts_summary.tsv)")
    return rows, unit_rows


def read_busco_json(path):
    """Return (lineage_name, completeness_float_or_None) from a BUSCO short_summary JSON."""
    with open(path) as fh:
        data = json.load(fh)
    results = data.get("results", {})
    m = re.search(r"C:\s*([\d.]+)%", results.get("one_line_summary", "") or "")
    completeness = float(m.group(1)) if m else results.get("Complete percentage")
    lineage = data.get("lineage_dataset", {}).get("name", "NA")
    return lineage, completeness


def collect_busco(cat_dir):
    """Read BUSCO JSONs from busco_lineage/ and busco_eukaryote/ under cat_dir.

    JSON names follow <species>_<taxonID>_<src>_Lbusco.json / _Ebusco.json where
    <src> is the annotation origin the parameters were trained on (e.g. 'reference',
    'lyric', 'isoquant').

    Returns:
      flat_rows  : list of dicts keyed by BUSCO_FLAT_COLS, one per (species, src)
      srcs_seen  : sorted list of unique src values found (e.g. ['isoquant', 'reference'])
    """
    rows = {}
    for path in sorted(glob.glob(os.path.join(cat_dir, "busco_lineage", "*_Lbusco.json"))):
        stem = os.path.basename(path)[: -len("_Lbusco.json")]
        sp, src = stem.rsplit("_", 1)
        lineage, completeness = read_busco_json(path)
        rows.setdefault((sp, src), {}).update(
            {"lineage_used": lineage, "lineage_completeness": completeness})
    for path in sorted(glob.glob(os.path.join(cat_dir, "busco_eukaryote", "*_Ebusco.json"))):
        stem = os.path.basename(path)[: -len("_Ebusco.json")]
        sp, src = stem.rsplit("_", 1)
        _, completeness = read_busco_json(path)
        rows.setdefault((sp, src), {})["eukaryote_completeness"] = completeness
    flat_rows = [
        {"species": sp, "src": src, **v}
        for (sp, src), v in sorted(rows.items())
    ]
    srcs_seen = sorted({src for (_, src) in rows})
    return flat_rows, srcs_seen


def write_busco_and_general(cat_dir, counts_rows):
    """Write busco_summary.tsv and general_summary.tsv for one evaluation flavour.

    busco_summary.tsv  : flat table, one row per (species, src) pair.
    general_summary.tsv: one row per species; all counts columns followed by BUSCO
                         columns pivoted by src (reference_lineage_completeness, ...).
                         lineage_used depends on the taxon, not the source, so the
                         first src seen for a species fills it. For flavours without
                         BUSCO (e.g. refined), only the counts columns are written.

    counts_rows is the per-species list of metric lists returned by write_counts.
    """
    flat_rows, srcs_seen = collect_busco(cat_dir)

    # busco_summary.tsv: flat, one row per (species, src)
    busco_path = os.path.join(cat_dir, "busco_summary.tsv")
    with open(busco_path, "w") as fh:
        fh.write("\t".join(BUSCO_FLAT_COLS) + "\n")
        for row in flat_rows:
            fh.write("\t".join(str(row.get(c, "NA")) for c in BUSCO_FLAT_COLS) + "\n")
    print(f"  busco  : {len(flat_rows)} entries -> busco_summary.tsv")

    # pivot BUSCO by src for general_summary
    busco_by_sp = {}
    for row in flat_rows:
        sp, src = row["species"], row["src"]
        d = busco_by_sp.setdefault(sp, {})
        d[f"{src}_lineage_completeness"] = row.get("lineage_completeness", "NA")
        d[f"{src}_eukaryote_completeness"] = row.get("eukaryote_completeness", "NA")
        # lineage_used depends on the taxon, not the source: first src seen fills it
        if "lineage_used" not in d:
            d["lineage_used"] = row.get("lineage_used", "NA")

    busco_gcols = (
        ["lineage_used"]
        + [f"{src}_{m}" for src in srcs_seen
           for m in ("lineage_completeness", "eukaryote_completeness")]
    ) if srcs_seen else []

    general_cols = METRIC_HEADER + busco_gcols
    counts_by_sp = {r[0]: r for r in counts_rows}
    all_species = sorted(set(counts_by_sp) | set(busco_by_sp))

    general_path = os.path.join(cat_dir, "general_summary.tsv")
    with open(general_path, "w") as fh:
        fh.write("\t".join(general_cols) + "\n")
        for sp in all_species:
            cnt = list(counts_by_sp.get(sp, [sp] + ["NA"] * (len(METRIC_HEADER) - 1)))
            bco = [str(busco_by_sp.get(sp, {}).get(c, "NA")) for c in busco_gcols]
            fh.write("\t".join(cnt + bco) + "\n")
    print(f"  general: {len(all_species)} species -> general_summary.tsv")


LOCAL_COLS = (["unit"] + BUSCO_FLAT_COLS[:1] + ["src"] + METRIC_HEADER[1:]
              + ["lineage_used", "lineage_completeness", "eukaryote_completeness"])


def write_local_summary(cat_dir, unit_rows):
    """Write general_summary_local.tsv: a flat table with one row per (species, source)
    unit, joining per-unit counts with that unit's BUSCO completeness.

    This table is for local inspection only (per-source detail); it is NOT part of the
    cross-project overall summary, which stays one row per species in general_summary.tsv.

    unit_rows is the list of dicts returned by write_counts.
    """
    flat_rows, _ = collect_busco(cat_dir)
    busco_by_unit = {(r["species"], r["src"]): r for r in flat_rows}

    path = os.path.join(cat_dir, "general_summary_local.tsv")
    with open(path, "w") as fh:
        fh.write("\t".join(LOCAL_COLS) + "\n")
        for u in unit_rows:
            b = busco_by_unit.get((u["sp"], u["src"]), {})
            row = ([u["unit"], u["sp"], u["src"]] + u["metrics"][1:]
                   + [str(b.get("lineage_used", "NA")),
                      str(b.get("lineage_completeness", "NA")),
                      str(b.get("eukaryote_completeness", "NA"))])
            fh.write("\t".join(row) + "\n")
    print(f"  local  : {len(unit_rows)} units -> general_summary_local.tsv")


def main():
    parser = argparse.ArgumentParser(
        description="Write gene/transcript model counts and genome size for each "
        "evaluation flavour (regular, merged) into results/summary/<flavour>/."
    )
    repo_root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
    default_results = os.path.join(repo_root, "results")
    default_species = os.path.join(repo_root, "training_data", "species")
    parser.add_argument("-r", "--results", default=os.path.normpath(default_results),
                        help="Path to the results/ directory (default: ../results relative to this script)")
    parser.add_argument("-g", "--species-root", default=os.path.normpath(default_species),
                        help="Directory holding per-species <sp>/CLEAN_*.fna genomes "
                        "(default: ../training_data/species relative to this script)")
    args = parser.parse_args()

    results_dir = args.results
    species_root = args.species_root
    summary_dir = os.path.join(results_dir, "summary")
    os.makedirs(summary_dir, exist_ok=True)

    counts = {}
    for name, cfg in CATEGORIES.items():
        counts[name] = write_counts(results_dir, summary_dir, name, cfg, species_root)

    print(f"\nDone. Counts written under {summary_dir}")

    for name in CATEGORIES:
        cat_dir = os.path.join(summary_dir, name)
        rows, unit_rows = counts[name]
        write_busco_and_general(cat_dir, rows)
        write_local_summary(cat_dir, unit_rows)

    for name in CATEGORIES:
        for kind in ("busco_lineage", "busco_eukaryote"):
            busco_dir = os.path.join(summary_dir, name, kind)
            run_busco_plot(os.path.join(busco_dir, "*.json"),
                           os.path.join(busco_dir, "busco_summary.png"))


if __name__ == "__main__":
    main()
