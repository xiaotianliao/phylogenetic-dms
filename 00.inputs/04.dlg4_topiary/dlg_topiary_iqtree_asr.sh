#!/bin/bash
#BSUB -G team354
#BSUB -n 20
#BSUB -M 200000
#BSUB -R "select[mem>200000] rusage[mem=200000]"
#BSUB -o /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/04.dlg4_topiary/00.log/dlg_topiary_asr_090424_%J.out
#BSUB -e /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/04.dlg4_topiary/00.log/dlg_topiary_asr_090424_%J.err
#BSUB -q normal

#Akaike Information Criterion:           Q.plant+R6
#Corrected Akaike Information Criterion: LG
#Bayesian Information Criterion:         Q.plant+R6
#Best-fit model: Q.plant+R6 chosen according to BIC

conda activate /software/hgi/envs/conda/team354/xl7/asr

iqtree2 -s /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/04.dlg4_topiary/seed_to_ali/07_alignment_trimmed.fasta -m Q.plant+R6 -asr -nt AUTO -pre /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/04.dlg4_topiary/02.out_asr/dlg_topiary_asr/dlg_topiary_asr



