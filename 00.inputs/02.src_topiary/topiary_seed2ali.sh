#!/bin/bash
#BSUB -G team354
#BSUB -n 10
#BSUB -M 300000
#BSUB -R "select[mem>300000] rusage[mem=300000]"
#BSUB -o /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/02.src_topiary/topiary_src_090224_%J.out
#BSUB -e /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/02.src_topiary/topiary_src_090224_%J.err
#BSUB -q hugemem

conda activate  /software/hgi/envs/conda/team354/xl7/topiary

cd /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/02.src_topiary

topiary-seed-to-alignment /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/02.src_topiary/seed-dataframe_src_tec_abl.csv --out_dir  /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/02.src_topiary/seed_to_ali  