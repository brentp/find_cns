import sys
import collections
from bx.intervals.intersection import IntervalTree, Interval


# brachy_v1   Bd2 4440720 4440745 Bradi2g05930    ATCCCCAAATTGTACTTGCACTAACA  sorghum_v1.4    3   2939646 2939689 Sb03g002930 TCTTGTTGGTAGTGCAAGTACGATTTGGTGGTTTGTTACTTTTA    rice_v6 1   5130702 5130727 Os01g09880  ATCCCCAAATTGTACTTGCACTAACA  c[at][at][ag]ttg[tg][at][ac][tg]tgca[ac][tg][at]ac

def parse_current_cnss(fn, org):
    """
    just grab rice_v6 from the triplets.
    add to interval tree.
    """
    tree = collections.defaultdict(IntervalTree)
    cnss = {}
    for i, sline in enumerate(open(fn)):
        if i == 0: continue # header
        line = sline.rstrip().split("\t")
        if not org in line: continue
        for i in range(0, 15, 6):
            if line[i] == org: break
        seqid, start, end, accn = line[i+1: i + 5]
        start, end = sorted(map(int, (start, end)))
        cnsid = (seqid, start, end, accn)
        iv = Interval(start, end, value=(accn, cnsid))

        tree[seqid].add_interval(iv)
        cnss[cnsid] = (seqid, start, end, accn)
    return tree, cnss



def main(cnsblast, tree, cnss, link, print_extra=False):
    seen = {}
    print "\t".join(["qseqid", "qstart", "qend", "qaccn", "sseqid", "sstart", "send", "link"])
    for sline in open(cnsblast):
        line = sline.split("\t")
        qseqid, sseqid = line[:2]
        qstart, qend, sstart, send = map(int, line[6:10])

        overlaps = tree[qseqid].find(qstart, qend)
        if len(overlaps) == 0:
            continue
        smid = (sstart + send) / 2
        for o in overlaps:
            accn, cnsid = o.value
            qmid = (o.start + o.end)/2
            line = map(str, [qseqid, o.start, o.end, accn, sseqid, sstart, send, link % locals()])
            line = "\t".join(line)
            if line in seen: continue
            try:
                del cnss[cnsid]
            except KeyError:
                pass
            seen[line] = 1
            print line
    if print_extra:
        for cnsid, (seqid, start, stop, accn) in cnss.iteritems():
            print "\t".join(map(str, (seqid, start, stop, accn)))


if __name__ == "__main__":

    current_cnss = sys.argv[1]
    current_org = sys.argv[2]
    tree, cnss = parse_current_cnss(current_cnss, current_org)
    cnsblast = sys.argv[3]

    current_dsgid = 9109
    newgsid = 10536

    link = "http://genomevolution.org/CoGe/GEvo.pl?prog=blastn;dsgid1=" + str(current_dsgid)
    link += ";chr1=%(qseqid)s;dsgid2=" + str(newgsid) + ";num_seqs=2;chr2=%(sseqid)s;"
    link += "x1=%(qmid)i;x2=%(smid)i;autogo=1"
    main(cnsblast, tree, cnss, link, print_extra=False)
