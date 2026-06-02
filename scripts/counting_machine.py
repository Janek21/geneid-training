import os, glob, argparse

#Counts gene/transcript models and creates a summary table
CATEGORIES = {
    "regular": {"pred_dir": "pred",         "pred_suffix": "_own.gff3"},
    "refined": {"pred_dir": "refined_pred", "pred_suffix": "_own.gff3"},
    "merged":  {"pred_dir": "merged",       "pred_suffix": "_own.gff"},
}


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


def write_counts(results_dir, summary_dir, name, cfg):
    """Write per-species gene/transcript model counts for one evaluation flavour into results/summary/<name>/: a counts/ folder with <sp>_gc.txt / <sp>_tc.txt integer files, plus a counts_summary.tsv table at the category root."""
    cat_out = os.path.join(summary_dir, name)
    counts_out = os.path.join(cat_out, "counts")
    os.makedirs(counts_out, exist_ok=True)

    pred_dir = os.path.join(results_dir, cfg["pred_dir"])
    suffix = cfg["pred_suffix"]
    rows = []
    for gff in sorted(glob.glob(os.path.join(pred_dir, "*" + suffix))):
        sp = os.path.basename(gff)[: -len(suffix)]  # <specie_name>_<taxid>
        genes, transcripts = count_models(gff)
        with open(os.path.join(counts_out, f"{sp}_gc.txt"), "w") as fh:
            fh.write(f"{genes}\n")
        with open(os.path.join(counts_out, f"{sp}_tc.txt"), "w") as fh:
            fh.write(f"{transcripts}\n")
        rows.append((sp, genes, transcripts))

    with open(os.path.join(cat_out, "counts_summary.tsv"), "w") as fh:
        fh.write("species\tgene_count\ttranscripts_count\n")
        for sp, genes, transcripts in rows:
            fh.write(f"{sp}\t{genes}\t{transcripts}\n")
    print(f"[{name}] counts: {len(rows)} species -> counts/ (+ counts_summary.tsv)")
    return rows


def main():
    parser = argparse.ArgumentParser(
        description="Write gene/transcript model counts for each evaluation flavour "
        "(regular, refined, merged) into results/summary/<flavour>/."
    )
    default_results = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
    parser.add_argument("-r", "--results", default=os.path.normpath(default_results),
                        help="Path to the results/ directory (default: ../results relative to this script)")
    args = parser.parse_args()

    results_dir = args.results
    summary_dir = os.path.join(results_dir, "summary")
    os.makedirs(summary_dir, exist_ok=True)

    for name, cfg in CATEGORIES.items():
        write_counts(results_dir, summary_dir, name, cfg)

    print(f"\nDone. Counts written under {summary_dir}")


if __name__ == "__main__":
    main()
