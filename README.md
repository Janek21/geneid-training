# geneid-training

Species-specific [geneid](https://github.com/guigolab/geneid) training, genome-wide prediction, and annotation merging for protist genomes. The workflow is driven by three Jupyter notebooks with SLURM array jobs for the heavy steps. Trained models and predictions can be evaluated against BUSCO and compared with the base annotation.

Engine repository - ab initio gene prediction trained on long-read evidence. Used alongside [LyRic_annotator](https://github.com/Janek21/LyRic_annotator) (unguided) and [isoquant_annotator](https://github.com/Janek21/isoquant_annotator) (guided) in a larger protist annotation pipeline.

## Overview

Each species/source unit (e.g. `Babesia_duncani_323732_lyric`) flows through three notebooks in order:

1. **`notebooks/1-GeneIDTrain.ipynb`** - cleans the assembly, extracts CDS from the base annotation (reference, LyRic, or IsoQuant), samples up to 1000 genes, and writes `job_commands/total_training.txt` (one geneid training block per unit).
2. **`notebooks/2-GeneIDPredict.ipynb`** - runs geneid with the trained parameters across each assembly to produce GFF3 predictions, evaluates them with BUSCO, and writes gene/transcript count metrics.
3. **`notebooks/3-AnnotationMerging.ipynb`** - merges each geneid prediction with its source CDS annotation using AGAT, evaluates the merged annotation with BUSCO, and writes merged-annotation count metrics.

The notebooks write command files to `job_commands/`; the SLURM wrappers in `scripts/memory/` dispatch these as arrays via `scripts/array_lib.sh`.

## Repository structure

```
geneid-training/
├── notebooks/
│   ├── 1-GeneIDTrain.ipynb       # prepare training data + write training command file
│   ├── 2-GeneIDPredict.ipynb     # run geneid + BUSCO on predictions
│   ├── 3-AnnotationMerging.ipynb # merge prediction + source annotation, re-evaluate
│   └── deprecated/               # earlier single-merged-annotation experiments
├── scripts/
│   ├── array_lib.sh              # shared SLURM array runner (sourced by memory/ wrappers)
│   ├── clean_ref.sh              # strips FASTA headers to bare sequence IDs
│   ├── get_CDS.sh                # extracts CDS features from a GFF/GFF3 (gz or plain)
│   ├── sample_CDS.sh             # randomly samples N genes from a CDS GFF
│   ├── missing_train.sh          # finds units in total_training.txt not yet trained
│   ├── CDS_machine.py            # CDS GFF → geneid-ready input (Python helpers)
│   ├── counting_machine.py       # gene/transcript model counting + summary TSV generation
│   ├── get_busco_db.py           # resolves BUSCO lineage per taxon (NCBI Entrez)
│   ├── get_genetic_code.py       # resolves NCBI translation table per taxon
│   ├── buscoPlot.py              # plots BUSCO completeness from result JSONs
│   └── memory/                   # SLURM wrappers (sbatch from job_commands/)
│       ├── trainer.sh            # array job for total_training.txt
│       ├── predictor.sh          # array job for total_preds.txt
│       ├── merger.sh             # array job for mergeAnn.txt
│       └── busco.sh              # array job for buscoScoring.txt or M_buscoScoring.txt
```

## Requirements

- SLURM cluster
- geneid installed and on `$PATH`
- conda env `buscomania`: AGAT, BUSCO, pigz
- Jupyter (for running the notebooks)
- NCBI Entrez access (email set inside the notebooks) for BUSCO lineage and genetic code lookups
- Reference data one level up from the repository(2 levels up from the notebooks):
  `../../data/species/<species_taxid>/GC[AF]*/` (genome FASTA + annotation GFF)
- Engine siblings checked out at the same level (`../../lyric_annotator/`, `../../isoquant_annotator/`) for their per-species output annotations

Optionally, clone the geneid reference parameter files for comparison:
```bash
git clone https://github.com/guigolab/geneid-parameter-files training_data/geneid-parameter-files
```

## Usage

### 1. Run the notebooks in order

Open and run `notebooks/1-GeneIDTrain.ipynb`, `notebooks/2-GeneIDPredict.ipynb`, and `notebooks/3-AnnotationMerging.ipynb` in sequence. Each notebook:
- Processes all species/source units defined at the top of the notebook
- Writes a command file to `job_commands/`

### 2. Dispatch the SLURM array jobs

Each SLURM wrapper in `scripts/memory/` is run from the `job_commands/` directory. When called without an array index it auto-counts the blocks and resubmits itself as an array:

```bash
cd job_commands

# train geneid parameters (one task per species/source unit)
sbatch trainer.sh

# predict genome-wide (one task per unit)
sbatch predictor.sh

# BUSCO on predictions
sbatch busco.sh buscoScoring.txt

# merge prediction with source CDS annotation (one task per unit)
sbatch merger.sh

# BUSCO on merged annotations
sbatch busco.sh M_buscoScoring.txt
```

If some training jobs failed, identify and resubmit only the missing units:
```bash
cd job_commands
bash ../scripts/missing_train.sh

#resubmit
sbatch trainer.sh missing_training.txt
```

### 3. Build the summary tables

Run the summary cells at the bottom of notebooks 2 and 3 (they call `counting_machine.write_counts`), or equivalently:

```python
# regular predictions
write_counts(results_dir, summary_dir, "regular", CATEGORIES["regular"])

# merged annotations
write_counts(results_dir, summary_dir, "merged", CATEGORIES["merged"])
```

## Output

```
training_data/species/<species>/
    CLEAN_*.fna             # header-stripped genome
    CDS_<unit>.gff          # extracted CDS features
    sample1000_CDS_<unit>.gff  # 1000-gene training sample

results/
    trainedParams/          # one .geneid.optimized.param file per unit
    pred/                   # geneid GFF3 predictions, one per unit
    merged/                 # AGAT-merged (prediction + source CDS) GFF, one per unit
    summary/
        regular/            # evaluation of geneid predictions
            busco_lineage/  # taxon-driven BUSCO JSONs
            busco_eukaryote/ # Eukaryota BUSCO JSONs
            counts/         # per-unit metrics TSVs
        merged/             # evaluation of merged annotations (same structure)
        counts_summary.tsv, busco_summary.tsv, general_summary.tsv  # aggregate tables

job_commands/
    total_training.txt      # geneid training commands (from notebook 1)
    total_preds.txt         # geneid prediction commands (from notebook 2)
    buscoScoring.txt        # BUSCO commands for predictions (from notebook 2)
    mergeAnn.txt            # AGAT merge commands (from notebook 3)
    M_buscoScoring.txt      # BUSCO commands for merged annotations (from notebook 3)
```
