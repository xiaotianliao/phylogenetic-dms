#!/bin/bash
#BSUB -G team354
#BSUB -n 20
#BSUB -M 400000
#BSUB -R "select[mem>400000] rusage[mem=400000]"
#BSUB -o /lustre/scratch126/gengen/projects/alpha-allostery-global/phylo-dms/02.iqtree_intermediates/00.log/src_homo_primates_asr_082724_%J.out
#BSUB -e /lustre/scratch126/gengen/projects/alpha-allostery-global/phylo-dms/02.iqtree_intermediates/00.log/src_homo_primates_asr_082724_%J.err
#BSUB -q hugemem

conda activate /software/hgi/envs/conda/team354/xl7/asr

cd /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/02.intermediates/03.foldmason_outs/01.all_structures

foldmason easy-msa /lustre/scratch126/gengen/projects/alpha-allostery-global/git-phylo-dms/02.intermediates/00.pdbs/00.all_human_primate/*.pdb result.fasta tmpFolder