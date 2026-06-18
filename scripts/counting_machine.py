import os, glob, argparse, sys

#Counts gene/transcript models and creates a summary table
_SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))
if _SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, _SCRIPTS_DIR)
from buscoPlot import run as _busco_plot


def run_busco_plot(glob_pattern, output_path):
    _busco_plot(glob_pattern, output_path)


CATEGORIES = {
    "regular": {"pred_dir": "pred",         "pred_suffix": "_own.gff3"},
    "refined": {"pred_dir": "refined_pred", "pred_suffix": "_own.gff3"},
    "merged":  {"pred_dir": "merged",       "pred_suffix": "_own.gff"},
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
    """Write per-species counts + derived metrics for one evaluation flavour into results/summary/<name>/: a counts/ folder with one consolidated <sp>_metrics.tsv per species, plus a counts_summary.tsv table (same columns) at the category root. Genome size is the total assembly length read from <species_root>/<sp>/CLEAN_*.fna; species_root defaults to <results_dir>/../training_data/species."""
    if species_root is None:
        species_root = os.path.join(os.path.dirname(os.path.abspath(results_dir)),
                                    "training_data", "species")
    cat_out = os.path.join(summary_dir, name)
    counts_out = os.path.join(cat_out, "counts")
    os.makedirs(counts_out, exist_ok=True)

    pred_dir = os.path.join(results_dir, cfg["pred_dir"])
    suffix = cfg["pred_suffix"]
    rows = []
    for gff in sorted(glob.glob(os.path.join(pred_dir, "*" + suffix))):
        sp = os.path.basename(gff)[: -len(suffix)]  # <specie_name>_<taxid>
        genes, transcripts = count_models(gff)
        gsize = genome_size_bp(species_root, sp)
        gd, td, ipg = derived_metrics(genes, transcripts, gsize)
        row = [sp, str(genes), str(transcripts),
               "NA" if gsize is None else str(gsize),
               "NA", "NA", gd, td, ipg, "NA"]
        with open(os.path.join(counts_out, f"{sp}_metrics.tsv"), "w") as fh:
            fh.write("\t".join(METRIC_HEADER) + "\n")
            fh.write("\t".join(row) + "\n")
        rows.append(row)

    with open(os.path.join(cat_out, "counts_summary.tsv"), "w") as fh:
        fh.write("\t".join(METRIC_HEADER) + "\n")
        for row in rows:
            fh.write("\t".join(row) + "\n")
    print(f"[{name}] metrics: {len(rows)} species -> counts/ (+ counts_summary.tsv)")
    return rows


def main():
    parser = argparse.ArgumentParser(
        description="Write gene/transcript model counts and genome size for each "
        "evaluation flavour (regular, refined, merged) into results/summary/<flavour>/."
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

    for name, cfg in CATEGORIES.items():
        write_counts(results_dir, summary_dir, name, cfg, species_root)

    print(f"\nDone. Counts written under {summary_dir}")

    # Two BUSCO runs per species (taxon lineage + eukaryota_odb12) land in
    # separate folders; plot each independently (mirrors the isoquant evaluation).
    for name in ("regular", "merged"):
        for kind in ("busco_lineage", "busco_eukaryote"):
            busco_dir = os.path.join(summary_dir, name, kind)
            run_busco_plot(os.path.join(busco_dir, "*.json"),
                           os.path.join(busco_dir, "busco_summary.png"))


if __name__ == "__main__":
    main()
