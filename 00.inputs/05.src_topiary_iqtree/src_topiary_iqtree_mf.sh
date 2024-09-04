#!/bin/bash
#BSUB -G team354
#BSUB -n 24
#BSUB -M 300000
#BSUB -R "select[mem>300000] rusage[mem=300000]"
#BSUB -o /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/05.src_topiary_iqtree/00.log/src_topiary_iqtree_mf_090324_%J.out
#BSUB -e /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/05.src_topiary_iqtree/00.log/src_topiary_iqtree_mf_090324_%J.err
#BSUB -q normal

conda activate /software/hgi/envs/conda/team354/xl7/asr

iqtree2 -s /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/02.src_topiary/seed_to_ali/08_alignment_cleaned_trimmed.fasta -m MF -cmax 8 -nt AUTO -pre /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/05.src_topiary_iqtree/01.out_mf/src_topiary_mf/src_topiary_mf

