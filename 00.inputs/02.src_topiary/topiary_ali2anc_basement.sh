#!/bin/bash
#BSUB -G team354
#BSUB -n 28  # Requesting 28 threads (cores)
#BSUB -M 300000  # Memory requirement (adjust as needed)
#BSUB -R "select[mem>300000] rusage[mem=300000] span[ptile=28]" 
#BSUB -o /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/02.src_topiary/ali2anc_src_basement_090324_%J.out
#BSUB -e /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/02.src_topiary/ali2anc_src_basement_090324_%J.err
#BSUB -q basement


conda activate  /software/hgi/envs/conda/team354/xl7/topiary

cd /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/02.src_topiary

topiary-alignment-to-ancestors /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/02.src_topiary/seed_to_ali/final_dataframe.csv --out_dir /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/02.src_topiary/ali_to_anc_basement --num_threads 28
