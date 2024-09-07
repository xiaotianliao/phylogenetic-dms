#!/bin/bash
#BSUB -G team354
#BSUB -n 20  # Requesting 10 threads (cores)
#BSUB -M 100000  # Memory requirement (adjust as needed)
#BSUB -R "select[mem>100000] rusage[mem=100000]" 
#BSUB -R "affinity[core(20)]"
#BSUB -o /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/02.src_topiary/ali2anc_src_090624_week_%J.out
#BSUB -e /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/02.src_topiary/ali2anc_src_090624_week_%J.err
#BSUB -q week

#topiary-fasta-into-dataframe 05_clean-aligned-dataframe.csv 10_aligned_cleaned_trimed_sh3_sh2_kd.fasta final_dataframe_10.csv

conda activate  /software/hgi/envs/conda/team354/xl7/topiary

cd /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/02.src_topiary

topiary-alignment-to-ancestors /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/02.src_topiary/seed_to_ali/final_dataframe_10.csv --out_dir /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/02.src_topiary/ali_to_anc_week --num_threads 20

