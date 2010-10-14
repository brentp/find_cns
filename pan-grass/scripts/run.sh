# run bsr.py then as_groups.py then glamify.py
Q=rice_v6
S=musa_v1.0

# sh scripts/quota.sh data/${Q} data/musa_v1.0 4:8
#python scripts/write_masked_genomic_fasta.py data/${Q}.bed data/${Q}.fasta > data/${Q}.masked.fasta
#python scripts/write_masked_genomic_fasta.py data/${S}.bed data/${S}.fasta > data/${S}.masked.fasta

#python /usr/local/src/bio-pipeline/blast_nearby/blast_nearby.py \
#                    --qbed data/${Q}.nolocaldups.bed \
#                    --sbed data/${S}.nolocaldups.bed \
#                    --anchors data/${Q}_${S}.pairs.txt \
#                    --dist 14000 \
#                    --cmd "bl2seq -p blastn -D 1 -E 2 -q -2 -r 1 -G 5 -W 7 -F T \
#                          -Y 812045000 -d 26195 -e 2.11 \
#                          -i %(query_fasta)s -j %(subject_fasta)s \
#                           | grep -v '#' | grep -v 'WARNING' | grep -v 'ERROR'" \
#                    data/${Q}.masked.fasta \
#                    data/${S}.masked.fasta > data/cnshits.blast

python scripts/merge_cnss.py data/bsr.with.motifs.triplets.csv rice_v6 data/cnshits.blast > data/rice_v6_musa_v1.0.pan-grass.csv
