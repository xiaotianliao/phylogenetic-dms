#!/bin/bash
#BSUB -G team354
#BSUB -n 20
#BSUB -M 300000
#BSUB -R "select[mem>300000] rusage[mem=300000]"
#BSUB -o /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/05.src_topiary_iqtree/00.log/src_topiary_asr_090324_%J.out
#BSUB -e /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/05.src_topiary_iqtree/00.log/src_topiary_asr_090324_%J.err
#BSUB -q normal

#Akaike Information Criterion:           Q.insect+I+R7
#Corrected Akaike Information Criterion: LG+G4
#Bayesian Information Criterion:         Q.insect+I+R7
#Best-fit model: Q.insect+I+R7 chosen according to BIC

conda activate /software/hgi/envs/conda/team354/xl7/asr

iqtree2 -s /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/02.src_topiary/seed_to_ali/08_alignment_cleaned_trimmed.fasta -m Q.insect+I+R7 -asr -nt AUTO -pre /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/05.src_topiary_iqtree/02.out_asr/src_topiary_asr/src_topiary_asr


