# write a genomic fasta file with all sequences covered by features
# in the specified Bed file masked to N.
from flatfeature import Bed
import sys
b = Bed(sys.argv[1], sys.argv[2])

for seqid, seq in b.mask_cds():
    print ">%s" % seqid
    print seq.tostring()
