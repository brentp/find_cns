from flatfeature import Bed
import sys
import os.path as op

bed = sys.argv[1]
genomic_fasta = sys.argv[2]

b = Bed(bed, genomic_fasta)
outfile = "%s.features%s" % op.splitext(genomic_fasta)

b.cds_fasta(outfile=outfile)
