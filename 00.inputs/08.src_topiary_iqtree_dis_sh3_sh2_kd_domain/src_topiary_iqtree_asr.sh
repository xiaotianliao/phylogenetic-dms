#!/bin/bash
#BSUB -G team354
#BSUB -n 10
#BSUB -M 100000
#BSUB -R "select[mem>100000] rusage[mem=100000]"
#BSUB -o /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/08.src_topiary_iqtree_dis_sh3_sh2_kd_domain/00.log/src_topiary_asr_090624_%J.out
#BSUB -e /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/08.src_topiary_iqtree_dis_sh3_sh2_kd_domain/00.log/src_topiary_asr_090624_%J.err
#BSUB -q yesterday

#Akaike Information Criterion:           Q.insect+F+I+R8
#Corrected Akaike Information Criterion: LG+G4
#Bayesian Information Criterion:         Q.insect+F+I+R8
#Best-fit model: Q.insect+F+I+R8 chosen according to BIC

conda activate /software/hgi/envs/conda/team354/xl7/asr

iqtree2 -s /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/02.src_topiary/seed_to_ali/09_alignment_cleaned_disordered_sh3_sh2_kd.fasta -m Q.insect+F+I+R8 -asr -nt AUTO -pre /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/08.src_topiary_iqtree_dis_sh3_sh2_kd_domain/02.out_asr/src_topiary_asr/src_topiary_asr


