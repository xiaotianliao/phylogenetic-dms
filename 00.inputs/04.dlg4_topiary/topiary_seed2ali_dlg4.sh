#!/bin/bash
#BSUB -G team354
#BSUB -n 10
#BSUB -M 100000
#BSUB -R "select[mem>100000] rusage[mem=100000]"
#BSUB -o /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/04.dlg4_topiary/topiary_dlg4_090324_%J.out
#BSUB -e /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/04.dlg4_topiary/topiary_dlg4_090324_%J.err
#BSUB -q yesterday

conda activate  /software/hgi/envs/conda/team354/xl7/topiary

cd /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/04.dlg4_topiary

topiary-seed-to-alignment /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/04.dlg4_topiary/seed-dataframe_dlg4.csv --out_dir  /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/00.inputs/04.dlg4_topiary/seed_to_ali  --local_recip_blast_db /lustre/scratch126/gengen/projects/alpha-allostery-global/phylo-dms/01.blastdb/dlg4_localblast_homo_danio/combined_protein.faa  