#!/bin/bash
#BSUB -G team354
#BSUB -n 24
#BSUB -M 100000
#BSUB -R "select[mem>100000] rusage[mem=100000]"
#BSUB -o /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/04.dlg4_topiary/00.log/dlg_topiary_iqtree_mf_090424_%J.out
#BSUB -e /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/04.dlg4_topiary/00.log/dlg_topiary_iqtree_mf_090424_%J.err
#BSUB -q normal

conda activate /software/hgi/envs/conda/team354/xl7/asr

iqtree2 -s /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/04.dlg4_topiary/seed_to_ali/07_alignment_trimmed.fasta -m MF -cmax 8 -nt AUTO -pre /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/04.dlg4_topiary/01.out_mf/dlg_topiary_mf/dlg_topiary_mf


