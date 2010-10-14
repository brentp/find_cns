from bx.intervals.intersection import Intersecter, Interval
import collections
import sys


def parse_datasheet(f, q_or_s, as_tree=False, sep=","):
    """
    parse a cns datasheet into an interval tree. q or s specifies wether
    to use the query or subject as the interval start/stop.
    """
    t = collections.defaultdict(Intersecter)
    fh = open(f)
    header = fh.readline().strip().split(sep)
    intable = ('qstart', 'qstop', 'sstart', 'sstop')
    for sline in fh:
        line = sline.strip().split(sep)
        d = dict(zip(header, line))
        if d['orthologous'] != 'T': continue

        for ib in intable: d[ib] = int(d[ib])
        if as_tree:
            t[d[q_or_s + 'chr']].add_interval(Interval(d[q_or_s + 'start'], d[q_or_s + 'stop'], d))
        else:
            yield d
    if as_tree:
        yield dict(t)
       
def main(fa, fb, org, akey="s", bkey="s"):
    tmpl = "%(org)s,%(seqid)s,%(accn)s,%(start)d,%(end)d"
    atree = parse_datasheet(fa, akey, as_tree=True).next()
    btree = parse_datasheet(fb, bkey, as_tree=True).next()
    out = open("%s.merged.csv" % org, "w")
    print >>sys.stderr, "writing to", out.name
    for b in parse_datasheet(fb, bkey, as_tree=False):
        accns = [b[bkey + 'accn']]
        t = btree[b[bkey + 'chr']]
        srng = set(range(b[bkey + 'start'], b[bkey + 'stop'] + 1))
        # need to change this in case it's not s
        overlaps = t.find(b[bkey + 'start'], b[bkey + 'stop'])
        if not overlaps: continue

        for intersected in (ii.value for ii in overlaps):
            this_srange = range(intersected[bkey + 'start'], intersected[bkey + 'stop'] + 1)
            srng = srng.union(this_srange)
            accns.append(intersected[bkey + 'accn'])

            #this_qrange = range(intersected['qstart'], intersected['qstop'] + 1)
            #qrng = qrng.union(this_qrange)

        srng = sorted(srng)
        print >> out, tmpl % dict(org=org, seqid=b[bkey + 'chr'], accn="|".join(sorted(set(accns))), start=srng[0], end=srng[-1])
if __name__ == "__main__":
    f1 = "/var/www/cns/v2.00/brachy_v1_sorghum_v1.4/brachy_v1_sorghum_v1.4/brachy_v1_sorghum_v1.4_2010-06-14_cns_datasheet.attrs.csv"
    f2 = "/var/www/cns/v2.00/rice_v6_sorghumv1.4/rice_v6_sorghum_v1.4/rice_v6_sorghum_v1.4_2010-06-17_cns_datasheet.attrs.csv"
    f3 = "/var/www/cns/v2.00/rice_v6_brachy_v1/rice_v6_brachy_v1/rice_v6_brachy_v1_2010-08-25_cns_datasheet.attrs.csv"

    main(f1, f2, "sorghum", "s", "s")
    main(f2, f3, "rice", "q", "q")
    main(f1, f3, "brachy", "q", "s")

