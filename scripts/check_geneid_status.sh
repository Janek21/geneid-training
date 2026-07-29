#!/bin/bash
# Classify every unit listed in units.txt by how far it got through geneid-training
# (notebooks 1 -> 2 -> 3).
# A "unit" is <species>_<source>, e.g. Babesia_duncani_323732_lyric
# (source is one of reference/lyric/isoquant, and never contains "_", so the species
# name is recovered as unit minus its last "_"-segment: unit.rsplit("_", 1)[0], same
# as the notebooks themselves do).
#
# units.txt is the full set of species_source units given as initial input to
# notebook 1 (one per line, plain "<species>_<source>", no quoting) -- supplied
# externally, not derived from any pipeline output. This matters because it's the
# only complete record of what should exist: total_training.txt (notebook 1's last
# cell) only lists units that already have a sampled CDS file by that point, so a
# unit that fails at or before CDS extraction never appears there at all. units.txt
# fixes that by being the pipeline's actual starting point rather than a downstream
# byproduct.
#
# Buckets (mutually exclusive, in pipeline order, mirroring notebook 1 -> 2 -> 3):
#   NO GENOME CLEANED   - species-level: clean_ref.sh never produced CLEAN_*.fna for this species
#   NO CDS EXTRACTED    - get_CDS.sh never produced (or produced empty, later auto-deleted) CDS_<unit>.gff
#   NOT SAMPLED         - sample_CDS.sh never produced sample1000_CDS_<unit>.gff
#   NOT TRAINED         - geneidTRAINer never produced/moved the optimized param file
#   TRAINED, NOT PRED   - trained param file exists, but no genome-wide GFF3 prediction yet
#   PRED, NOT EVAL      - prediction exists, but BUSCO (lineage + eukaryote) hasn't fully run
#   EVAL, NOT MERGED    - prediction BUSCO-evaluated, but not yet AGAT-merged with source annotation
#   MERGED, NOT EVAL    - merged annotation exists, but BUSCO on it hasn't fully run
#   MERGED + EVALUATED  - merged annotation exists and BUSCO-evaluated (fully complete)
#
# Stage artifacts (confirmed exact filenames against Janek21/geneid-training, not
# globbed guesses; verified against notebooks 1/2/3 on 2026-07-28):
#   clean genome  training_data/species/<species>/CLEAN_*.fna                          (per species)
#   CDS extracted training_data/species/<species>/CDS_<unit>.gff                       (per unit)
#   sampled       training_data/species/<species>/sample1000_CDS_<unit>.gff            (per unit)
#   trained       results/trainedParams/<unit>.geneid.optimized.param
#   predicted     results/pred/<unit>.gff3
#   pred eval     results/summary/regular/busco_lineage/<unit>_Lbusco.json AND
#                 results/summary/regular/busco_eukaryote/<unit>_Ebusco.json   (both required)
#   merged        results/merged_pred/<unit>.gff
#   merge eval    results/summary/merged/busco_lineage/<unit>_Lbusco.json AND
#                 results/summary/merged/busco_eukaryote/<unit>_Ebusco.json   (both required)
#
# NOTE: merged output lives in results/merged_pred/, not results/merged/. Notebook 3
# does `mkdir -p ../results/merged_pred` and `--out ../results/merged_pred/<unit>.gff`.
# (scripts/counting_machine.py's CATEGORIES["merged"] still says "merged" -- looks like
# a separate, pre-existing bug there, independent of this checker.)
#
# Each unit's BUSCO step always launches two jobs (lineage + eukaryote) together, so
# "eval" here requires BOTH json files -- a unit with only one is genuinely partway
# through, not done.
#
# Usage (run from the repo root, or pass paths explicitly):
#   bash check_geneid_status.sh [units.txt] [base_dir]
# Defaults:
#   units.txt = units.txt
#   base_dir  = .   (where training_data/, results/ live)

shopt -s nullglob

units_file="${1:-units.txt}"
base_dir="${2:-.}"
training_data_dir="$base_dir/training_data/species"
results_dir="$base_dir/results"
summary_dir="$base_dir/results/summary"

sample_n=1000  # matches n=1000 in notebook 1's sampling + training-command cells

if [ ! -s "$units_file" ]; then
	echo "units.txt not found or empty: $units_file" >&2
	exit 1
fi

# species = unit minus its last "_"-delimited segment (the source, e.g. lyric/isoquant/reference)
species_of() {
	printf '%s\n' "${1%_*}"
}

#--- stage detection ---

has_clean_genome() {  # $1 = species
	local f=("$training_data_dir/$1"/CLEAN_*.fna)
	[ "${#f[@]}" -gt 0 ] && [ -s "${f[0]}" ]
}

has_cds() {  # $1 = species, $2 = unit
	local f="$training_data_dir/$1/CDS_$2.gff"
	[ -s "$f" ]
}

has_sampled() {  # $1 = species, $2 = unit
	local f="$training_data_dir/$1/sample${sample_n}_CDS_$2.gff"
	[ -s "$f" ]
}

has_trained() {  # $1 = unit
	local f="$results_dir/trainedParams/$1.geneid.optimized.param"
	[ -s "$f" ]
}

has_pred() {  # $1 = unit
	local f="$results_dir/pred/$1.gff3"
	[ -s "$f" ]
}

has_pred_eval() {  # $1 = unit -- requires BOTH lineage and eukaryote BUSCO output
	local l="$summary_dir/regular/busco_lineage/$1_Lbusco.json"
	local e="$summary_dir/regular/busco_eukaryote/$1_Ebusco.json"
	[ -s "$l" ] && [ -s "$e" ]
}

has_merged() {  # $1 = unit
	local f="$results_dir/merged_pred/$1.gff"
	[ -s "$f" ]
}

has_merged_eval() {  # $1 = unit -- requires BOTH lineage and eukaryote BUSCO output
	local l="$summary_dir/merged/busco_lineage/$1_Lbusco.json"
	local e="$summary_dir/merged/busco_eukaryote/$1_Ebusco.json"
	[ -s "$l" ] && [ -s "$e" ]
}

#--- classify ---

no_clean=()
no_cds=()
not_sampled=()
not_trained=()
trained_only=()
pred_not_eval=()
eval_not_merged=()
merged_only=()
merged_eval=()

while IFS= read -r unit || [ -n "$unit" ]; do
	unit="${unit%$'\r'}"                 # strip stray CR
	[ -z "${unit// }" ] && continue      # skip blank lines
	[[ "$unit" == \#* ]] && continue     # skip comment lines

	sp=$(species_of "$unit")

	if ! has_clean_genome "$sp"; then
		no_clean+=("$unit")
	elif ! has_cds "$sp" "$unit"; then
		no_cds+=("$unit")
	elif ! has_sampled "$sp" "$unit"; then
		not_sampled+=("$unit")
	elif ! has_trained "$unit"; then
		not_trained+=("$unit")
	elif ! has_pred "$unit"; then
		trained_only+=("$unit")
	elif ! has_pred_eval "$unit"; then
		pred_not_eval+=("$unit")
	elif ! has_merged "$unit"; then
		eval_not_merged+=("$unit")
	elif ! has_merged_eval "$unit"; then
		merged_only+=("$unit")
	else
		merged_eval+=("$unit")
	fi
done < "$units_file"

#--- report ---

print_group() {
	local title="$1"; shift
	printf '\n== %s (%d) ==\n' "$title" "$#"
	if [ "$#" -eq 0 ]; then
		echo "  (none)"
	else
		printf '  %s\n' "$@"
	fi
}

total=$(( ${#no_clean[@]} + ${#no_cds[@]} + ${#not_sampled[@]} + ${#not_trained[@]} \
	+ ${#trained_only[@]} + ${#pred_not_eval[@]} + ${#eval_not_merged[@]} \
	+ ${#merged_only[@]} + ${#merged_eval[@]} ))
echo "Units file:    $units_file"
echo "Base dir:      $base_dir"
echo "Total units:   $total"

print_group "NO GENOME CLEANED (species-level, blocks everything)" "${no_clean[@]}"
print_group "NO CDS EXTRACTED"                                      "${no_cds[@]}"
print_group "CDS EXTRACTED but NOT SAMPLED"                         "${not_sampled[@]}"
print_group "SAMPLED but NOT TRAINED"                                "${not_trained[@]}"
print_group "TRAINED but NOT PREDICTED"                              "${trained_only[@]}"
print_group "PREDICTED but NOT EVALUATED (BUSCO incomplete)"         "${pred_not_eval[@]}"
print_group "EVALUATED but NOT MERGED"                               "${eval_not_merged[@]}"
print_group "MERGED but NOT EVALUATED (BUSCO incomplete)"            "${merged_only[@]}"
print_group "MERGED + EVALUATED (complete)"                          "${merged_eval[@]}"
