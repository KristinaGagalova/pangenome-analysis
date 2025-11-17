# index variants
bgzip bnapus16.distances.tsv.edges.weights.txt.community.0.fa.gz.a65af12.11fba48.e6193c9.smooth.final.Darmor.vcf
tabix bnapus16.distances.tsv.edges.weights.txt.community.0.fa.gz.a65af12.11fba48.e6193c9.smooth.final.Darmor.vcf.gz

# create index based on genome
vg autoindex --workflow giraffe -r BnapusDarmor-bzh_chromosomes.fasta -v bnapus16.distances.tsv.edges.weights.txt.community.0.fa.gz.a65af12.11fba48.e6193c9.smooth.final.Darmor.vcf.gz -p bnapus16.distances.tsv.edges.weights.txt.community.0.fa.gz.a65af12.11fba48.e6193c9.smooth.final.Darmor
