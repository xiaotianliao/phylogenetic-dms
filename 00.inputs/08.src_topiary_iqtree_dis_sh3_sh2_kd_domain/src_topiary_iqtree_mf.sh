#!/bin/bash
#BSUB -G team354
#BSUB -n 10
#BSUB -M 100000
#BSUB -R "select[mem>100000] rusage[mem=100000]"
#BSUB -o /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/08.src_topiary_iqtree_dis_sh3_sh2_kd_domain/00.log/src_topiary_iqtree_mf_090624_%J.out
#BSUB -e /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/08.src_topiary_iqtree_dis_sh3_sh2_kd_domain/00.log/src_topiary_iqtree_mf_090624_%J.err
#BSUB -q yesterday

conda activate /software/hgi/envs/conda/team354/xl7/asr

iqtree2 -s /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/02.src_topiary/seed_to_ali/09_alignment_cleaned_disordered_sh3_sh2_kd.fasta -m MF -cmax 8 -nt AUTO -pre /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/08.src_topiary_iqtree_dis_sh3_sh2_kd_domain/01.out_mf/src_topiary_mf/src_topiary_mf

