import os
import pathlib
import glob
import subprocess

base_dir=pathlib.Path("../data/species")
tr_sp=[d.name for d in base_dir.iterdir() if d.is_dir()]
n=1000
#tr_sp=["Mytilus_galloprovincialis"]
os.makedirs("../job_commands", exist_ok=True)
with open("../job_commands/SAR_training.txt", "w") as out:

    for sp in tr_sp:

        # Use glob to find files instead of !realpath
        sampleN_list=list(base_dir.glob(f"{sp}/sample{n}*.gff"))
        clean_ref_list=list(base_dir.glob(f"{sp}/CLEAN*.fna"))
        
        if not sampleN_list or not clean_ref_list:
            print(f"Files not found for {sp}, skipping.")
            continue
            
        sampleN_CDS=sampleN_list[0].resolve()
        clean_ref=clean_ref_list[0].resolve()
        result_folder=pathlib.Path(f"../results/{sp}/")
        
        print(sampleN_CDS)
        print(clean_ref)
        print(sp)
        print(result_folder)

        #training command using singularity container
        geneidtrainer_command=f"singularity run \
                            ../data/geneidtrainerdocker.sif \
                            /scripts_geneid/geneidTRAINer4docker.pl \
                            -species {sp} \
                            -gff {sampleN_CDS} \
                            -fastas {clean_ref} \
                            -results {result_folder} \
                            -reduced no"
        geneidtrainer_command=geneidtrainer_command.replace("                             ", " ")
        print(geneidtrainer_command)
        out.write(geneidtrainer_command + '\n')
