#!/usr/bin/bash

module load singularity/4.1.0-slurm

MASH=/scratch/pawsey1142/kgagalova/pangenomes/partitioning/mash:2.3--he348c14_1
GENOMES=/scratch/pawsey1142/kgagalova/pangenomes/genomes/select_scaffolds/subset/all_genomes16.fasta.gz

singularity run $MASH mash dist $GENOMES $GENOMES -s 10000 -i > bnapus16.distances.tsv
python3 pggb/scripts/mash2net.py -m bnapus16.distances.tsv
python3.6 pggb/scripts/net2communities.py \
    -e bnapus16.distances.tsv.edges.list.txt \
    -w bnapus16.distances.tsv.edges.weights.txt \
    -n bnapus16.distances.tsv.vertices.id2name.txt

seq 0 19 | while read i; do
    chromosomes=$(cat bnapus16.distances.tsv.edges.weights.txt.community.$i.txt | cut -f 3 -d '#' | sort | uniq | tr '\n' ' ');
    echo "community $i --> $chromosomes";
done

ls bnapus16.distances.tsv.edges.weights.txt.community.*.txt | parallel -I{} bash extract_communities.sh ../genomes/select_scaffolds/subset/all_genomes16.fasta.gz {}
