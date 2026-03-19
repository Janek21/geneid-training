#gets all species which are not yet trained and reates an executable file to train them

#get all species which have a parameters file already
ls ../results/trainedParams/*.geneid.optimized.param | xargs -n 1 basename | sed 's/\.geneid\.optimized\.param//' > ../results/solved_species.txt
#get all lines from species that dont have a parameters file yet(species that are yet to be trained)
grep -Fvf ../results/solved_species.txt ../job_commands/total_training.txt > ../job_commands/missing_training.txt

#bash missing_training.txt

