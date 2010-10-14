org=musa_v1.0
other=rice_v6

#wget -O data/musa_v1.0.fasta http://genomevolution.org/CoGe/data/genomic_sequence/0/0/10/10536/10536.faa
# on syteny: perl export_to_bed.pl -dsg 10536 > musa_v1.0.bed

python scripts/write_cds_fasta.py data/${org}.bed data/${org}.fasta
python scripts/write_cds_fasta.py data/${other}.bed data/${other}.fasta


