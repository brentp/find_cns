import sys
import os.path as op
sys.path.insert(0, op.abspath(op.dirname(__file__)))
from grouper import Grouper
from bsr import parse_datasheet
from pyfasta import Fasta

from bx.intervals.intersection import Intersecter, Interval
import collections

f1 = "brachy.merged.csv"
f2 = "sorghum.merged.csv"
f3 = "rice.merged.csv"

fa1 = Fasta('/home/brentp/work/cnspipeline/data/rice_v6_brachy_v1/brachy_v1.fasta')
fa2 = Fasta('/home/brentp/work/cnspipeline/data/rice_v6_sorghum_v1.4/sorghum_v1.4.fasta')
fa3 = Fasta('/home/brentp/work/cnspipeline/data/rice_v6_sorghum_v1.4/rice_v6.fasta')

fastas = {
    "brachy_v1": fa1,
    "sorghum_v1.4": fa2,
    "rice_v6": fa3
}


cns12 = "/var/www/cns/v2.00/brachy_v1_sorghum_v1.4/brachy_v1_sorghum_v1.4/brachy_v1_sorghum_v1.4_2010-06-14_cns_datasheet.attrs.csv"
cns32 = "/var/www/cns/v2.00/rice_v6_sorghumv1.4/rice_v6_sorghum_v1.4/rice_v6_sorghum_v1.4_2010-06-17_cns_datasheet.attrs.csv"
cns31 = "/var/www/cns/v2.00/rice_v6_brachy_v1/rice_v6_brachy_v1/rice_v6_brachy_v1_2010-08-25_cns_datasheet.attrs.csv"

# s, s
# q, q
# q, s

# main(cns12, cns23, "sorghum", "s", "s")
# main(cns23, cns21, "rice", "q", "q")
# main(cns12, cns21, "brachy", "q", "s")

def gen_cns(cns_file, qorg, sorg):
    header = None
    non = 0
    for sline in open(cns_file):
        if sline[0] == "#": continue
        line = sline.rstrip().split(",")
        if header is None:
            header = line
            continue
        cns = dict(zip(header, line))
        if cns['orthologous'] != 'T':
            non += 1
            continue
        sstart, sstop = sorted((int(cns['sstart']), int(cns['sstop'])))
        yield (cns['qchr'], int(cns['qstart']), int(cns['qstop']), qorg), \
              (cns['schr'], sstart, sstop, sorg)

    print >>sys.stderr, cns_file, "non:", non
def main():
    tree_brachy = merged_tree(f1)
    tree_sorghum = merged_tree(f2)
    tree_rice = merged_tree(f3)


    iter12 = gen_cns(cns12, "brachy_v1", "sorghum_v1.4")
    iter32 = gen_cns(cns32, "rice_v6", "sorghum_v1.4")
    iter31 = gen_cns(cns31, "rice_v6", "brachy_v1")

    groups = Grouper()

    """
    cns12 = "/var/www/cns/v2.00/brachy_v1_sorghum_v1.4/brachy_v1_sorghum_v1.4/brachy_v1_sorghum_v1.4_2010-06-14_cns_datasheet.attrs.csv"
    cns32 = "/var/www/cns/v2.00/rice_v6_sorghumv1.4/rice_v6_sorghum_v1.4/rice_v6_sorghum_v1.4_2010-06-17_cns_datasheet.attrs.csv"
    cns31 = "/var/www/cns/v2.00/rice_v6_brachy_v1/rice_v6_brachy_v1/rice_v6_brachy_v1_2010-08-25_cns_datasheet.attrs.csv"
    """

    # brachy-sorghum
    for qcns, scns in iter12:
        joiner(qcns, scns, tree_brachy, tree_sorghum, groups)

    # rice-sorghum
    for qcns, scns in iter32:
        joiner(qcns, scns, tree_rice, tree_sorghum, groups)

    # rice-brachy
    for qcns, scns in iter31:
        joiner(qcns, scns, tree_rice, tree_brachy, groups)

    n = 0
    for group in groups:
        if len(group) < 3: continue
        orgs = set(g[0] for g in group)
        if len(orgs) < 3: continue
        #print group
        #raw_input("\n... press any key ...\n")
        for org, seqid, start, end, accns in group:
            print (">" + org), seqid, start, end, accns
            print fastas[org][seqid][start - 1: end]
        n += 1
        print
    print >>sys.stderr, "groups of 3+:", n


def joiner(qcns, scns, qtree, stree, groups):
    qms = qtree[qcns[0]].find(qcns[1], qcns[2])
    if not qms: return
    sms = stree[scns[0]].find(scns[1], scns[2])
    if not sms: return

    for qm in qms:
        # org, seqid, start, end
        #print
        #print qcns[0], qm
        qhashable = (qcns[-1], qcns[0], qm.start, qm.end, qm.value)
        for sm in sms:
            #print "\t", scns[0], sm
            shashable = (scns[-1], scns[0], sm.start, sm.end, sm.value)
            #print qhashable, shashable
            #raw_input("\n... any key...\n")
            groups.join(qhashable, shashable)


def merged_tree(merged_f):
    trees = collections.defaultdict(Intersecter)
    for sline in open(merged_f):
        if sline[0] == "#": continue
        line = sline.rstrip().split(",")
        seqid, accns = line[1:3]
        start, end = map(int, line[3:])
        trees[seqid].add_interval(Interval(start, end, accns))

    return trees


if __name__ == "__main__":
    main()
